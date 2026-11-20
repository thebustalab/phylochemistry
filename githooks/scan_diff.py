#!/usr/bin/env python3
"""
scan_diff.py — the content gate behind githooks/pre-push.

    scan_diff.py --push LOCAL_SHA REMOTE_SHA REMOTE_NAME [--identity]

is what the hook runs, once per pushed ref. It checks EVERYTHING the push would send that is
not already on that remote — every commit in the range, not just the tip — and exits 1 with
findings on stderr if anything matches (0 when clean):

  * content: each commit's added lines (`git log -p`, so a secret added in one commit and
    deleted in the next is still caught). Hits are reported as `path:line @commit`, and flagged
    when the file is gitignored but still tracked;
  * binary blobs (and files a .gitattributes `-diff` hides from the patch), scanned raw;
  * commit messages, annotated-tag messages, and the paths of every new file;
  * with --identity (the public_tools repos): the identity grep, every author / committer /
    tagger email must be a GitHub noreply, and every image blob's metadata is parsed and must
    hold only allowlisted harmless fields (pixel size, resolution, orientation, colour profile,
    the comment "Screenshot") — GPS, make/model/serial, names, paths, edit history and any field
    the parser does not recognise all block, named but never printed.

Fails CLOSED: if the range cannot be computed (remote ref not fetched), the push is blocked
with a reason. The offending literal is printed partly redacted, so a genuine hit does not
paste the real secret into the terminal / a CI log.

TWO TIERS.  PATTERNS (always on) is the credential / private-path core.  IDENTITY_PATTERNS
(--identity) is public_tools/AGENTS.md's pre-publish grep
(`/home/bustalab|/Users/bust0037|d.umn.edu|umn.edu|@gmail`) plus the x500, made automatic. It
matches the BARE `umn.edu` domain, so a link to a umn.edu page fires too: a visible false
positive you bypass with SKIP_SECRET_SCAN=1 is the accepted price of the checklist being
enforced rather than remembered. It is opt-in because the academic WEBSITE repo legitimately
publishes a contact address, and a gate that fires on every push gets disabled. This file stays
byte-identical in every repo that carries it — the scope lives in the caller (the hook), not in
a forked copy.  Regression test: githooks/test_scan_diff.py.

Not checked: compressed content (zip/gzip/Office/PDF streams) and document or video metadata,
blobs over BLOB_SCAN_LIMIT, and anything inside the scanner's own three files (SELF_FILES).
`git push --no-verify` skips every hook, this one included.
"""
import re
import subprocess
import sys

PATTERNS = [
    (r"sk-[A-Za-z0-9_-]{20,}", "OpenAI/Anthropic-style secret key (sk-...)"),
    (r"ghp_[A-Za-z0-9]{30,}", "GitHub personal access token (ghp_...)"),
    (r"github_pat_[A-Za-z0-9_]{20,}", "GitHub fine-grained PAT (github_pat_...)"),
    (r"AKIA[0-9A-Z]{16}", "AWS access key ID (AKIA...)"),
    (r"AIza[0-9A-Za-z_\-]{30,}", "Google API key (AIza...)"),
    (r"xox[baprs]-[A-Za-z0-9\-]{10,}", "Slack token (xox...)"),
    (r"-----BEGIN [A-Z ]*PRIVATE KEY-----", "PEM private key block"),
    (r"export\s+\w*_API_KEY\w*\s*=\s*\S", "shell export of an *_API_KEY value"),
    (r"export\s+\w*_TOKEN\w*\s*=\s*\S", "shell export of an *_TOKEN value"),
    (r"\bACCESS_TOKEN\w*\s*=\s*[\"']", 'ACCESS_TOKEN = "..." literal'),
    (r"(?i)\b\w*API_?KEY\w*\s*=\s*[\"']", 'API_KEY = "..." literal'),
    # `.env` / config style: a secret-named key set to an unquoted or quoted mixed letter+digit value
    (r"(?i)\b\w*(?:api_?key|secret|token|passw(?:or)?d)\w*[\"']?\s*[:=]\s*[\"']?"
     r"(?=[^\s\"']*\d)(?=[^\s\"']*[A-Za-z])[A-Za-z0-9_\-~./+]{16,}(?![A-Za-z0-9_(])",
     "secret-named key set to a token-like value"),
    (r"\b\d{1,6}~[A-Za-z0-9]{40,}", "Canvas API token (NNNN~...)"),
    (r"(?i)\bBearer\s+[A-Za-z0-9._~+/-]{16,}", "HTTP Bearer token"),
    (r"[\"'][A-Fa-f0-9]{32,}[\"']", "long hex literal (32+ chars), token-shaped"),
    (r"[\"'][A-Za-z0-9+/]{40,}={0,2}[\"']", "long base64-shaped literal (40+ chars)"),
    # paths and hosts: case-insensitive, no trailing word boundary (/home/bustalab_old counts)
    (r"(?i)/home/bustalab", "private path (/home/bustalab)"),
    (r"(?i)/Users/bust0037", "private path (/Users/bust0037)"),
    (r"(?i)(?<![A-Za-z0-9])bustalab-desktop", "internal hostname (host2)"),
    (r"(?i)(?<![A-Za-z0-9])system76-pc", "internal hostname (host1)"),
    (r"(?i)(?<![A-Za-z0-9])spark-e5b5", "internal hostname (dgx)"),
    # bustalab.d.umn.edu is deliberately NOT here (Lucas, 2026-10-16): it is the lab server's
    # public name, already served in phylochemistry.R and analyzeGCMSdata6/7.R, which print the
    # Shiny "Connect at" URL from it. The --identity tier still catches it via d.umn.edu.
    (r"(?<![0-9.])131\.212\.57\.\d{1,3}", "internal lab IP (131.212.57.0/24)"),
]
# Random bytes in a binary blob match these two often enough to be noise; skip them there.
HEURISTIC = {"long hex literal (32+ chars), token-shaped", "long base64-shaped literal (40+ chars)"}
COMPILED = [(re.compile(p), label) for p, label in PATTERNS]

# The four assignment-shaped rules above match `NAME = <anything>`, so they also match a
# DOCUMENTED placeholder — `export CANVAS_TOKEN="…"` in canvas-assignments' README was a
# standing false positive that would have blocked a fresh-branch push. Only these four are
# placeholder-suppressible; a `sk-`/`ghp_`/PEM prefix match never is, because those shapes
# cannot be a placeholder.
PLACEHOLDER_SUPPRESSIBLE = {
    "shell export of an *_API_KEY value",
    "shell export of an *_TOKEN value",
    'ACCESS_TOKEN = "..." literal',
    'API_KEY = "..." literal',
    "secret-named key set to a token-like value",
}
_VALUE_RE = re.compile(r"""=\s*["'`]?([^"'`\n]*)""")
_PLACEHOLDER_PREFIX = ("REPLACE", "YOUR", "INSERT", "PASTE", "XXX", "CHANGEME", "TODO", "EXAMPLE")


def is_placeholder_value(tail: str) -> bool:
    """True when the value assigned in `tail` is an obvious stand-in, not a real secret."""
    m = _VALUE_RE.search(tail)
    if not m:
        return False
    v = m.group(1).strip()
    if not v or v in ("…", "...", "<...>"):
        return True
    if v[0] in "<${":
        return True
    if v.upper().startswith(_PLACEHOLDER_PREFIX):
        return True
    if re.fullmatch(r"[a-z_ ]+", v):             # `your_key`, `api_key`: words, not a secret
        return True
    return set(v.upper()) <= {"X", "*", ".", "-", "_", "…"}

# Opt-in tier (`--identity`), for the public TOOL repos. See the module docstring.
UMN_EMAIL = "email address (umn.edu)"
IDENTITY_PATTERNS = [
    (r"[A-Za-z0-9._%+-]*@[A-Za-z0-9.-]*umn\.edu", UMN_EMAIL),
    (r"umn\.edu", "UMN domain (checklist: bare umn.edu; a plain link is a false positive)"),
    (r"[A-Za-z0-9._%+-]*@gmail", "personal email address (@gmail)"),
    (r"bust0037", "Lucas's x500 (bust0037)"),
    (r"(?<![A-Za-z0-9])lucasbusta\d*@", "Lucas's personal email local-part"),
]
COMPILED_IDENTITY = [(re.compile(p, re.IGNORECASE), label) for p, label in IDENTITY_PATTERNS]


def redact(literal: str) -> str:
    if len(literal) <= 10:
        return literal[0] + "..." + literal[-1] if len(literal) > 2 else "***"
    return literal[:4] + "…redacted…" + literal[-4:]


def match_line(line: str, identity: bool, binary: bool = False):
    """Yield (label, redacted literal) for every rule `line` trips."""
    hit_email = False
    for regex, label in COMPILED + (COMPILED_IDENTITY if identity else []):
        if binary and label in HEURISTIC:
            continue
        if label in PLACEHOLDER_SUPPRESSIBLE:
            # flag if ANY match on the line is real (`A=YOUR_X; B=realsecret` used to pass)
            m = next((m for m in regex.finditer(line) if not is_placeholder_value(line[m.start():])), None)
        else:
            m = regex.search(line)
        if not m:
            continue
        if label == UMN_EMAIL:
            hit_email = True
        elif label.startswith("UMN domain") and hit_email:
            continue                             # the email hit already covers it
        yield label, redact(m.group(0))


# The scanner's own source IS a list of the patterns it looks for, so a diff that edits it
# trips every rule. That is not hypothetical — it blocked the very first commit of a freshly
# scaffolded tool. So exactly these three repo-root files are skipped, nothing else: a skip
# that matched any `githooks/` path let a secret through in `src/githooks/` (adversarial check,
# 2026-10-08). Accepted blind spot: a secret inside these three files.
SELF_FILES = {"githooks/scan_diff.py", "githooks/test_scan_diff.py", "githooks/pre-push"}

_HUNK_RE = re.compile(r"^@@ -\d+(?:,\d+)? \+(\d+)(?:,\d+)? @@")
_COMMIT_RE = re.compile(r"^\x00commit ([0-9a-f]+)$")
_BINARY_RE = re.compile(r"^Binary files (.+) and (.+) differ$")


def _unquote_path(path: str) -> str:
    # git C-quotes paths with unusual characters ("a b/c\303\251.txt"); keep it readable.
    if len(path) >= 2 and path[0] == path[-1] == '"':
        try:
            return path[1:-1].encode("latin-1").decode("unicode_escape").encode("latin-1").decode("utf-8")
        except (UnicodeError, ValueError):
            return path[1:-1]
    return path


def _strip_prefix(path: str) -> str:
    path = _unquote_path(path.strip())
    return path[2:] if path.startswith(("a/", "b/")) else path


def parse_patch(text: str):
    """Walk `git diff`/`git log -p` output.

    Returns (added, binary_paths): added is a list of (commit, path, new-file line, text);
    binary_paths are files whose content the patch did not show ("Binary files ... differ").
    Splits on '\\n' ONLY — str.splitlines() also breaks on form feeds and U+2028, which threw
    line numbers off and let a `+++ ` hidden after a U+2028 pose as a file header.
    """
    added, binary_paths = [], set()
    commit, path, lineno, in_header = "", "", 0, False
    for raw in text.split("\n"):
        m = _COMMIT_RE.match(raw)
        if m:
            commit, in_header = m.group(1), False
            continue
        if raw.startswith("diff --git "):
            in_header, path = True, ""
            continue
        if in_header:
            # file headers only count BEFORE the first hunk: an added line reading `++ x`
            # shows as `+++ x` inside a hunk, and used to switch the path to a skipped one
            if raw.startswith("+++ "):
                path = _strip_prefix(raw[4:])
            elif raw.startswith("--- ") and not path:
                path = _strip_prefix(raw[4:])
            elif _BINARY_RE.match(raw):
                b = _BINARY_RE.match(raw).group(2)
                binary_paths.add(_strip_prefix(b if b != "/dev/null" else _BINARY_RE.match(raw).group(1)))
            hunk = _HUNK_RE.match(raw)
            if hunk:
                in_header, lineno = False, int(hunk.group(1))
            continue
        hunk = _HUNK_RE.match(raw)
        if hunk:
            lineno = int(hunk.group(1))
            continue
        if raw.startswith("+"):
            if path not in SELF_FILES:
                added.append((commit, path, lineno, raw[1:]))
            lineno += 1
        elif raw.startswith(" "):
            lineno += 1
    return added, binary_paths


def scan(diff_text: str, identity: bool = False):
    """Return a list of (path, line, label, redacted literal) for a unified diff."""
    added, _ = parse_patch(diff_text)
    return [(p, n, label, sample) for _c, p, n, line in added for label, sample in match_line(line, identity)]


# ------------------------------------------------------------------ authors
NOREPLY_RE = re.compile(r"(@users\.noreply\.github\.com|^noreply@github\.com)$", re.IGNORECASE)


def scan_authors(log_text: str):
    """Flag author/committer emails (tab-separated `sha author committer`) that are not noreply.

    The content scan reads files; a commit's author field is separate, permanent and
    machine-harvestable once pushed (public_tools/AGENTS.md pre-publish checklist).
    """
    findings = []
    for raw in log_text.split("\n"):
        parts = raw.split("\t")
        if len(parts) < 3:
            continue
        sha, author, committer = parts[0], parts[1].strip(), parts[2].strip()
        for role, email in (("author", author), ("committer", committer)):
            if not NOREPLY_RE.search(email):
                findings.append((sha, role, redact(email) if email else "(empty)"))
    return findings


# ------------------------------------------------------------------ images
# The gate exists to stop PRIVATE information leaking through image metadata — location, device,
# serial numbers, owner/author names, host or user names, file paths, edit history. It is not a
# metadata-purity rule (Lucas, 2026-10-08): a macOS screenshot's pixel size, resolution, colour
# profile and the comment "Screenshot" are harmless and must not block a push. So every metadata
# block is PARSED and each field checked against an allowlist of harmless fields; anything else —
# including any field this parser does not recognise — blocks and is named. Unknown is not
# harmless. Findings name the field, never its value.
#
# Threat model: metadata a camera, phone, OS or editor writes without being asked. Deliberate
# steganography in the pixel data is out of scope.

# PNG chunks that carry pixels or colour, with the byte length each may have (None = any length
# within the bound). iDOT is Apple's decoding hint in every macOS screenshot (28 bytes of offsets).
PNG_PIXEL_CHUNKS = {
    b"IHDR": (13, 13), b"PLTE": (0, 768), b"IDAT": None, b"IEND": (0, 0), b"tRNS": (0, 256),
    b"cHRM": (32, 32), b"gAMA": (4, 4), b"sBIT": (1, 4), b"sRGB": (1, 1), b"bKGD": (1, 6),
    b"hIST": (0, 512), b"pHYs": (9, 9), b"acTL": (8, 8), b"fcTL": (26, 26), b"fdAT": None,
    b"iDOT": (28, 28), b"cICP": (4, 4), b"cLLI": (8, 8), b"mDCV": (24, 24),
}
# Metadata signatures searched for anywhere in an image. A hit that falls outside a block the
# parser has vetted means metadata is hiding where the walk did not look.
METADATA_SIGNATURES = [(b"Exif\x00\x00", "EXIF block"), (b"<x:xmpmeta", "XMP packet"),
                       (b"<?xpacket", "XMP packet"), (b"Photoshop 3.0\x00", "IPTC/Photoshop block")]

# ---- free text that is allowed in a few named fields
_SAFE_PUNCT = set(" .,:()_+-\u200e\u200f")     # the last two are bidi marks in Arabic/Hebrew names


def _safe_text(s: str) -> bool:
    """Short product/version-style text: letters (any script), digits and ` .,()_+-` only — so no
    path separator, `@`, `~` or quote — no possessive, and nothing the content rules flag.
    Used for Software / CreatorTool / the profile description."""
    import unicodedata
    s = s.strip()
    if not 0 < len(s) <= 64 or not all(c.isalnum() or c in _SAFE_PUNCT or
                                       unicodedata.category(c).startswith("M") for c in s):
        return False
    if re.search(r"(?i)\b\w+'s\b|\busers?\b|\bhome\b", s):
        return False
    # a hostname (`janes-mbp.local`, `lab.example.edu`) or an attribution (`Edited by Jane`)
    if re.search(r"(?i)\.(?:local|localdomain|lan|home|internal|corp|com|edu|org|gov|ac)\b|\bby\s+\w", s):
        return False
    return not any(True for _ in match_line(s, identity=True))


# ---- EXIF (TIFF-structured IFDs), as found in JPEG APP1, PNG eXIf and WebP EXIF
EXIF_NAMES = {
    0x0100: "ImageWidth", 0x0101: "ImageLength", 0x0102: "BitsPerSample", 0x0103: "Compression",
    0x0106: "PhotometricInterpretation", 0x010E: "ImageDescription", 0x010F: "Make",
    0x0110: "Model", 0x0111: "StripOffsets", 0x0112: "Orientation", 0x0115: "SamplesPerPixel",
    0x011A: "XResolution", 0x011B: "YResolution", 0x0128: "ResolutionUnit", 0x0131: "Software",
    0x0132: "DateTime", 0x013B: "Artist", 0x013C: "HostComputer", 0x0201: "JPEGInterchangeFormat",
    0x0202: "JPEGInterchangeFormatLength", 0x0213: "YCbCrPositioning", 0x8298: "Copyright",
    0x8769: "ExifIFD", 0x8825: "GPSInfo (GPS IFD)", 0x9000: "ExifVersion",
    0x9003: "DateTimeOriginal", 0x9004: "DateTimeDigitized", 0x9010: "OffsetTime",
    0x9011: "OffsetTimeOriginal", 0x9012: "OffsetTimeDigitized", 0x9101: "ComponentsConfiguration",
    0x927C: "MakerNote", 0x9286: "UserComment", 0x9290: "SubSecTime", 0x9291: "SubSecTimeOriginal",
    0x9292: "SubSecTimeDigitized", 0x9C9B: "XPTitle", 0x9C9C: "XPComment", 0x9C9D: "XPAuthor",
    0x9C9E: "XPKeywords", 0x9C9F: "XPSubject", 0xA000: "FlashpixVersion", 0xA001: "ColorSpace",
    0xA002: "PixelXDimension", 0xA003: "PixelYDimension", 0xA005: "InteropIFD",
    0xA420: "ImageUniqueID", 0xA430: "CameraOwnerName", 0xA431: "BodySerialNumber",
    0xA432: "LensSpecification", 0xA433: "LensMake", 0xA434: "LensModel",
    0xA435: "LensSerialNumber", 0xC4A5: "PrintIM", 0xEA1C: "Padding",
}
_TIFF_TYPE_SIZE = {1: 1, 2: 1, 3: 2, 4: 4, 5: 8, 6: 1, 7: 1, 8: 2, 9: 4, 10: 8, 11: 4, 12: 8}
_SHORT, _LONG, _RATIONAL, _UNDEF, _ASCII = {3}, {3, 4}, {5}, {7}, {2}
SCREENSHOT = "Screenshot"


def _exif_value_ok(tag, typ, count, raw, text):
    """Allowlist for IFD0 / Exif-IFD fields: (tag, type, count, value) must look harmless."""
    def one_number(types):
        return typ in types and count == 1
    if tag in (0x011A, 0x011B):                                   # X/YResolution
        return one_number(_RATIONAL)
    if tag == 0x0128:                                             # ResolutionUnit
        return one_number(_SHORT)
    if tag in (0x0112, 0x0213, 0xA001):                           # Orientation, YCbCrPositioning, ColorSpace
        return one_number(_SHORT)
    if tag in (0x0100, 0x0101, 0xA002, 0xA003):                   # image / pixel dimensions
        return one_number(_LONG)
    if tag in (0x9000, 0xA000):                                   # Exif / Flashpix version: 4 digits
        return typ in _UNDEF and count == 4 and raw[:4].isdigit()
    if tag == 0x9101:                                             # ComponentsConfiguration
        return typ in _UNDEF and count == 4 and all(b <= 6 for b in raw[:4])
    if tag == 0x9286:                                             # UserComment
        return typ in _UNDEF and (raw.rstrip(b"\x00 ") in (b"ASCII\x00\x00\x00" + SCREENSHOT.encode(),
                                                            b"\x00" * 8 + SCREENSHOT.encode()))
    if tag == 0x010E:                                             # ImageDescription
        return typ in _ASCII and text == SCREENSHOT
    if tag == 0x0131:                                             # Software: product/version only
        return typ in _ASCII and _safe_text(text)
    return False


def exif_fields(tiff: bytes):
    """Return the names of every EXIF field in `tiff` (a TIFF-structured block) that is not on
    the harmless allowlist. GPS, thumbnails (IFD1 — may show the uncropped original), MakerNote
    and anything unrecognised are always named."""
    bad = []
    if tiff[:4] == b"II*\x00":
        order = "little"
    elif tiff[:4] == b"MM\x00*":
        order = "big"
    else:
        return ["EXIF block (not TIFF-structured, cannot be parsed)"]

    def u(off, n):
        if off < 0 or off + n > len(tiff):
            raise ValueError("offset out of range")
        return int.from_bytes(tiff[off:off + n], order)

    seen = set()

    def walk(off, where, allowed_sub):
        if off in seen:
            bad.append(f"EXIF {where} (IFD loop)")
            return
        seen.add(off)
        n = u(off, 2)
        if n > 500:
            raise ValueError("implausible IFD")
        for k in range(n):
            e = off + 2 + 12 * k
            tag, typ, count = u(e, 2), u(e + 2, 2), u(e + 4, 4)
            name = EXIF_NAMES.get(tag, f"tag 0x{tag:04X}")
            size = _TIFF_TYPE_SIZE.get(typ, 0) * count
            if not size:
                bad.append(f"EXIF {name}")
                continue
            voff = e + 8 if size <= 4 else u(e + 8, 4)
            u(voff, size) if size <= 8 else u(voff + size - 1, 1)     # bounds check
            raw = tiff[voff:voff + size]
            if tag in allowed_sub:                                     # a sub-IFD pointer
                if typ not in _LONG | {13} or count != 1:
                    bad.append(f"EXIF {name}")
                else:
                    walk(u(voff, 4), allowed_sub[tag], {} if tag == 0xA005 else {0xA005: "InteropIFD"})
                continue
            if where == "InteropIFD":
                ok = tag in (0x0001, 0x0002) and size <= 4
            else:
                text = raw.split(b"\x00", 1)[0].decode("latin-1") if typ in _ASCII else ""
                ok = _exif_value_ok(tag, typ, count, raw, text)
                if ok and typ in _ASCII and raw[len(text) + 1:].strip(b"\x00"):
                    ok = False                       # text hidden after the terminating NUL
            if not ok:
                bad.append(f"EXIF {name}")
        nxt = u(off + 2 + 12 * n, 4)
        if nxt and where == "IFD0":
            bad.append("EXIF IFD1 (thumbnail — may show the uncropped original)")
        elif nxt:
            bad.append(f"EXIF {where} (chained IFD)")

    try:
        walk(u(4, 4), "IFD0", {0x8769: "ExifIFD"})
    except (ValueError, IndexError):
        bad.append("EXIF block (malformed, cannot be parsed)")
    return list(dict.fromkeys(bad))


# ---- XMP
_XMP_NS = {
    "http://ns.adobe.com/exif/1.0/": "exif", "http://ns.adobe.com/tiff/1.0/": "tiff",
    "http://ns.adobe.com/xap/1.0/": "xmp", "http://ns.adobe.com/photoshop/1.0/": "photoshop",
    "http://purl.org/dc/elements/1.1/": "dc", "http://ns.adobe.com/xap/1.0/mm/": "xmpMM",
    "http://ns.adobe.com/xap/1.0/sType/ResourceEvent#": "stEvt",
    "http://ns.adobe.com/xap/1.0/rights/": "xmpRights", "http://ns.adobe.com/exif/1.0/aux/": "aux",
    "http://cipa.jp/exif/1.0/": "exifEX", "http://iptc.org/std/Iptc4xmpCore/1.0/xmlns/": "Iptc4xmpCore",
    "http://ns.adobe.com/camera-raw-settings/1.0/": "crs", "http://ns.adobe.com/pdf/1.3/": "pdf",
    "http://ns.apple.com/adjustment-settings/1.0/": "AppleAdjust", "http://ns.google.com/photos/1.0/camera/": "GCamera",
}
_RDF = "http://www.w3.org/1999/02/22-rdf-syntax-ns#"
_XML = "http://www.w3.org/XML/1998/namespace"
_INT = re.compile(r"\d{1,6}")
_RATIO = re.compile(r"\d{1,7}(?:/\d{1,7})?(?:\.\d{1,4})?")
# property -> value check. Everything not here is named and blocks.
XMP_ALLOWED = {
    ("exif", "PixelXDimension"): _INT.fullmatch, ("exif", "PixelYDimension"): _INT.fullmatch,
    ("exif", "ColorSpace"): _INT.fullmatch, ("exif", "UserComment"): lambda v: v == SCREENSHOT,
    ("tiff", "ImageWidth"): _INT.fullmatch, ("tiff", "ImageLength"): _INT.fullmatch,
    ("tiff", "Orientation"): _INT.fullmatch, ("tiff", "ResolutionUnit"): _INT.fullmatch,
    ("tiff", "XResolution"): _RATIO.fullmatch, ("tiff", "YResolution"): _RATIO.fullmatch,
    ("photoshop", "ColorMode"): _INT.fullmatch, ("photoshop", "ICCProfile"): _safe_text,
    ("dc", "description"): lambda v: v == SCREENSHOT,
    ("dc", "format"): re.compile(r"image/[a-z0-9.+-]{1,20}").fullmatch,
    ("xmp", "CreatorTool"): _safe_text,
}


# toolkit name + up to four version-ish tokens, each containing a digit ("9.1-c002 79.dba3da3b5,
# 2023/12/15-10:42:37") — so no room for a word like a name
_XMPTK_RE = re.compile(r"(?:Adobe XMP Core|XMP Core|XMP toolkit|Image::ExifTool|Exiv2)"
                       r"(?:[ ,]+[0-9A-Za-z.\-/:+]*\d[0-9A-Za-z.\-/:+]*){1,4}[ ,]*")
_LANG_RE = re.compile(r"[A-Za-z]{1,8}(?:-[A-Za-z0-9]{1,8}){0,3}")



def xmp_fields(packet: bytes):
    """Names of every XMP property in `packet` that is not on the harmless allowlist."""
    import xml.etree.ElementTree as ET
    text = packet.decode("utf-8", "replace")
    if "<!DOCTYPE" in text or "<!ENTITY" in text:
        return ["XMP DOCTYPE/ENTITY (refused)"]
    # an xpacket wrapper is two processing instructions around the x:xmpmeta element
    if "<!--" in text:
        return ["XMP comment"]
    if re.search(r"<\?(?!xpacket\b)", text):
        return ["XMP processing instruction"]
    for uri in re.findall(r"""xmlns:?[\w.-]*\s*=\s*["']([^"']*)["']""", text):
        if uri not in _XMP_NS and uri not in (_RDF, "adobe:ns:meta/"):
            return ["XMP (unrecognised namespace declared)"]
    # the outermost element: first opening tag to LAST closing tag (find/rfind, not a lazy regex —
    # that was quadratic on a chunk of unclosed tags). Two packets inside it fail to parse.
    for open_, close in (("<x:xmpmeta", "</x:xmpmeta>"), ("<rdf:RDF", "</rdf:RDF>")):
        a, b = text.find(open_), text.rfind(close)
        if a != -1 and b > a:
            break
    else:
        return ["XMP packet (cannot be parsed)"]
    m_start, m_end = a, b + len(close)
    # outside the one packet: only the <?xpacket?> wrapper, whitespace and padding
    outside = re.sub(r"<\?xpacket[^>]{0,200}\?>", "", text[:m_start] + text[m_end:])
    if outside.strip(" \t\r\n\x00"):
        return ["XMP data outside the packet (a second packet or stray text)"]
    try:
        root = ET.fromstring(text[m_start:m_end])
    except ET.ParseError:
        return ["XMP packet (malformed, cannot be parsed)"]
    bad = []

    def qname(tag):
        ns, _, local = tag[1:].partition("}") if tag.startswith("{") else ("", "", tag)
        return _XMP_NS.get(ns), local

    def label(prefix, local):
        return f"XMP {prefix}:{local}" if prefix else f"XMP (unrecognised namespace):{local}"

    def simple_value(el):
        """The text of a simple property, or of an rdf:Alt holding exactly one x-default li."""
        kids = list(el)
        if not kids:
            return (el.text or "").strip()
        if len(kids) == 1 and kids[0].tag == f"{{{_RDF}}}Alt":
            lis = list(kids[0])
            if len(lis) == 1 and lis[0].tag == f"{{{_RDF}}}li" and not list(lis[0]) and \
                    set(lis[0].attrib) <= {f"{{{_XML}}}lang"}:
                return (lis[0].text or "").strip()
        return None

    def check(prefix, local, value):
        rule = XMP_ALLOWED.get((prefix, local))
        if not (rule and value is not None and rule(value)):
            bad.append(label(prefix, local))

    if root.tag == f"{{{_RDF}}}RDF":
        rdf = root
    else:
        kids = list(root)
        if len(kids) != 1 or kids[0].tag != f"{{{_RDF}}}RDF" or (root.text or "").strip() \
                or (kids[0].tail or "").strip():
            return ["XMP packet (not exactly one rdf:RDF)"]
        rdf = kids[0]
    if (rdf.text or "").strip() or set(rdf.attrib):
        bad.append("XMP rdf:RDF (stray text or attributes)")
    for desc in rdf:
        if (desc.text or "").strip() or (desc.tail or "").strip():
            bad.append("XMP rdf:Description (stray text)")
        if desc.tag != f"{{{_RDF}}}Description":
            bad.append(label(*qname(desc.tag)))
            continue
        for k, v in desc.attrib.items():
            if k == f"{{{_RDF}}}about":
                if v.strip():
                    bad.append("XMP rdf:about")
            else:
                check(*qname(k), v.strip())
        for prop in desc:
            prefix, local = qname(prop.tag)
            if set(prop.attrib) - {f"{{{_XML}}}lang"} or (prop.tail or "").strip() or \
                    not all(_LANG_RE.fullmatch(e.get(f"{{{_XML}}}lang", "en")) for e in prop.iter()):
                bad.append(label(prefix, local) + " (attributes or stray text)")
                continue
            check(prefix, local, simple_value(prop))
    for k, v in root.attrib.items():                         # x:xmptk ("XMP Core 6.0.0") only
        if not k.endswith("}xmptk") or not _XMPTK_RE.fullmatch(v.strip()):
            bad.append(f"XMP {qname(k)[1]} attribute")
    return list(dict.fromkeys(bad))


# ---- ICC colour profiles (PNG iCCP, JPEG APP2, WebP ICCP)
# Tags that are colourimetry, curves or profile bookkeeping: numbers, no provenance.
ICC_NUMERIC_TAGS = {
    b"wtpt", b"bkpt", b"rXYZ", b"gXYZ", b"bXYZ", b"rTRC", b"gTRC", b"bTRC", b"kTRC", b"chad",
    b"chrm", b"lumi", b"meas", b"tech", b"view", b"A2B0", b"A2B1", b"A2B2", b"B2A0", b"B2A1",
    b"B2A2", b"gamt", b"pre0", b"pre1", b"pre2", b"cicp", b"rig0", b"ciis", b"clro", b"clrt",
    b"clot", b"D2B0", b"D2B1", b"D2B2", b"B2D0", b"B2D1", b"B2D2",
    b"aarg", b"aagg", b"aabg", b"vcgt", b"vcgp", b"ndin", b"arts",    # Apple / ArgyllCMS curves
}
ICC_TEXT_TAGS = {b"desc", b"dscm"}                 # the profile description (Apple: multi-language)
# The profile COPYRIGHT is free text a calibration tool may fill with a person's name, so only
# known vendor/standard strings pass.
_ICC_COPYRIGHT_RE = re.compile(
    r"\s*(?:Copyright|\(c\)|©)?\s*(?:\(c\)|©)?\s*(?:\d{4}(?:-\d{4})?,?\s*)?(?:Apple Inc\.?|"
    r"Apple Computer,? Inc\.?|Adobe Systems,? Inc(?:orporated|\.)?|International Color Consortium|"
    r"Microsoft Corporation|Hewlett-Packard Company|Google,? (?:Inc\.|LLC)|Canon Inc\.?|"
    r"Samsung(?: Electronics(?: Co\., Ltd\.)?)?|Sony Corporation|Eastman Kodak Company|"
    r"X-Rite,? Inc\.?|FUJIFILM Corporation|Nikon Corporation)(?:[\s.,]*\d{4}(?:-\d{4})?)?"
    r"(?:[\s.,]*all rights reserved)?[\s.,]*|"
    r"No copyright, use freely|"
    r"Public Domain|CC0|FB", re.I)
_SKIA_DESC = re.compile(r"Google/Skia/[0-9A-Fa-f]{32}")       # Chrome's generated profiles
# Apple's 'mmod' (display make + model + serial + manufacture date). Make and model are the same
# information as the profile description, which is allowed; a display SERIAL is not. 0xFD626D62 is
# not a per-unit serial: it is the value every Apple-silicon built-in panel reports (seen on
# unrelated MacBooks, e.g. displayplacer's `s4251086178`), and it is what both published covers
# carry. Any other non-zero serial, or a manufacture date, blocks.
APPLE_SHARED_PANEL_SERIALS = {0, 0xFD626D62}


def _icc_text(data: bytes):
    """Decode an ICC text-ish tag (desc / mluc / text); None if the type is not understood."""
    t = data[:4]

    def cstr(b):
        text, _, rest = b.partition(b"\x00")
        return None if rest.strip(b"\x00") else text.decode("latin-1")
    if t == b"text":
        return [cstr(data[8:])]
    if t == b"desc":                         # ASCII part, then a UTF-16 part, then a Mac script part
        n = int.from_bytes(data[8:12], "big")
        u = 12 + n
        un = int.from_bytes(data[u + 4:u + 8], "big")
        uni = data[u + 8:u + 8 + 2 * un].decode("utf-16-be", "replace").rstrip("\x00")
        sc = u + 8 + 2 * un
        out = [cstr(data[12:u]), uni, cstr(data[sc + 3:sc + 70])]
        if len(data) > sc + 70 and data[sc + 70:].strip(b"\x00"):
            out.append(None)
        return out
    if t == b"mluc":
        n, rec = int.from_bytes(data[8:12], "big"), int.from_bytes(data[12:16], "big")
        out = []
        for k in range(min(n, 200)):
            r = 16 + rec * k
            ln, off = int.from_bytes(data[r + 4:r + 8], "big"), int.from_bytes(data[r + 8:r + 12], "big")
            out.append(data[off:off + ln].decode("utf-16-be", "replace"))
        return out
    return None


def icc_fields(prof: bytes):
    """Names of every ICC-profile tag not on the harmless allowlist."""
    if len(prof) < 132 or prof[36:40] != b"acsp":
        return ["ICC profile (not a valid profile, cannot be parsed)"]
    bad = []
    n = int.from_bytes(prof[128:132], "big")
    if n > 200 or 132 + 12 * n > len(prof):
        return ["ICC profile (malformed tag table)"]
    for k in range(n):
        sig = prof[132 + 12 * k:136 + 12 * k]
        off = int.from_bytes(prof[136 + 12 * k:140 + 12 * k], "big")
        size = int.from_bytes(prof[140 + 12 * k:144 + 12 * k], "big")
        data = prof[off:off + size]
        name = "ICC " + (sig.decode("latin-1") if sig.isalnum() or sig.strip().isalnum() else sig.hex())
        if off + size > len(prof):
            bad.append(name + " (out of range)")
        elif sig in ICC_NUMERIC_TAGS:
            continue
        elif sig in ICC_TEXT_TAGS:
            texts = _icc_text(data)
            if texts is None or None in texts or not all(_safe_text(t) or _SKIA_DESC.fullmatch(t)
                                                          for t in texts if t.strip()):
                bad.append(name + " (profile description)")
        elif sig == b"cprt":
            texts = _icc_text(data)
            if texts is None or None in texts or not all(_ICC_COPYRIGHT_RE.fullmatch(t) for t in texts):
                bad.append(name + " (copyright text)")
        elif sig == b"mmod" and data[:4] == b"mmod" and len(data) >= 24:
            serial, date = int.from_bytes(data[16:20], "big"), int.from_bytes(data[20:24], "big")
            if serial not in APPLE_SHARED_PANEL_SERIALS or date or data[24:].strip(b"\x00"):
                bad.append("ICC mmod (display serial number / manufacture date)")
        else:
            bad.append(name + (" (device make/model text)" if sig in (b"dmnd", b"dmdd") else ""))
    return list(dict.fromkeys(bad))


# ---- IPTC: no IPTC field is harmless (by-line, city, caption, keywords...), so it always blocks;
# the parse exists only to name the datasets.
IPTC_NAMES = {5: "ObjectName", 25: "Keywords", 55: "DateCreated", 60: "TimeCreated", 80: "By-line",
              85: "By-lineTitle", 90: "City", 92: "Sub-location", 95: "Province-State",
              100: "Country-Code", 101: "Country", 105: "Headline", 110: "Credit", 115: "Source",
              116: "CopyrightNotice", 118: "Contact", 120: "Caption", 122: "Writer-Editor"}


def iptc_fields(irb: bytes):
    """Name the Photoshop resources / IPTC datasets in an APP13 'Photoshop 3.0' block."""
    bad, i = [], irb.find(b"8BIM")
    while 0 <= i and i + 12 <= len(irb) and irb[i:i + 4] == b"8BIM":
        rid = int.from_bytes(irb[i + 4:i + 6], "big")
        nlen = irb[i + 6]
        j = i + 6 + 1 + nlen + ((1 + nlen) & 1)
        size = int.from_bytes(irb[j:j + 4], "big")
        body = irb[j + 4:j + 4 + size]
        if rid == 0x0404:
            p = 0
            while p + 5 <= len(body) and body[p] == 0x1C:
                rec, ds, ln = body[p + 1], body[p + 2], int.from_bytes(body[p + 3:p + 5], "big")
                bad.append(f"IPTC {rec}:{ds} {IPTC_NAMES.get(ds, '') if rec == 2 else ''}".rstrip())
                p += 5 + ln
        else:
            bad.append(f"Photoshop resource 0x{rid:04X}")
        i = j + 4 + size + (size & 1)
    return list(dict.fromkeys(bad)) or ["IPTC/Photoshop block"]


HEIF_BRANDS = {b"heic", b"heix", b"hevc", b"hevx", b"heim", b"heis", b"hevm", b"hevs", b"heif",
               b"mif1", b"mif2", b"mif3", b"msf1", b"avif", b"avis", b"avio", b"avci", b"avcs",
               b"crx ", b"jxl ", b"jpx ", b"jp2 "}


def image_kind(data: bytes):
    if data.startswith(b"\x89PNG\r\n\x1a\n"):
        return "png"
    if data.startswith(b"\xff\xd8"):
        return "jpeg"
    if data[:4] == b"RIFF" and data[8:12] == b"WEBP":
        return "webp"
    if data[:6] in (b"GIF87a", b"GIF89a"):
        return "gif"
    if data[:4] in (b"II*\x00", b"MM\x00*", b"II+\x00", b"MM\x00+"):   # TIFF, BigTIFF
        return "tiff"
    if data[4:8] == b"ftyp":
        box = data[8:min(int.from_bytes(data[:4], "big"), 256)]
        brands = {box[k:k + 4] for k in range(0, len(box) - 3, 4)}       # major, version, compatible
        if brands & HEIF_BRANDS:
            return "heif"
    if data[:4] == b"8BPS" or data[:9] == b"gimp xcf " or data[:2] == b"\xff\x0a" \
            or data[:12] in (b"\x00\x00\x00\x0cjP  \r\n\x87\n", b"\x00\x00\x00\x0cJXL \r\n\x87\n"):
        return "unreadable"                                             # PSD, XCF, JPEG XL, JPEG 2000
    return None


_PNG_TEXT_KEYWORDS = {"Title", "Author", "Description", "Copyright", "Creation Time", "Software",
                      "Disclaimer", "Warning", "Source", "Comment", "XML:com.adobe.xmp",
                      "Raw profile type exif", "Raw profile type iptc", "Raw profile type xmp",
                      "date:create", "date:modify", "date:timestamp", "parameters", "prompt"}


def _png_text_keyword(ctype, body):
    kw = body.split(b"\x00", 1)[0].decode("latin-1")
    return f"PNG {ctype} '{kw}'" if kw in _PNG_TEXT_KEYWORDS else f"PNG {ctype} (unrecognised keyword)"


XMP_SIGS = [b"<x:xmpmeta", b"<?xpacket", b"<?xpacket"]     # one packet, begin + end wrapper


def image_metadata(data: bytes, path: str = ""):
    """Return the private/unknown metadata fields in an image (empty = clean, or not an image).

    Detection is by CONTENT, not extension, so a TIFF renamed .png is still caught.
    """
    import zlib
    kind = image_kind(data)
    found, vetted = [], []                     # vetted: (start, end) byte spans already parsed

    def vet(start, end, fields, sigs=()):
        """Record a parsed block's findings and the signatures its parse accounts for (each
        allowed at most once — a second packet hiding in the same span is still caught)."""
        found.extend(fields)
        vetted.append((start, end, list(sigs)))

    if kind in ("tiff", "heif", "unreadable"):
        return ["TIFF/HEIC/AVIF/PSD/XCF/JPEG XL/JPEG 2000 cannot be checked here; "
                "export to PNG/JPEG and strip"]
    if kind == "png":
        i = 8
        while i + 8 <= len(data):
            n = int.from_bytes(data[i:i + 4], "big")
            ctype = data[i + 4:i + 8]
            body = data[i + 8:i + 8 + n]
            cname = ctype.decode("latin-1") if ctype.isalpha() else ctype.hex()
            if i + 12 + n > len(data):
                found.append(f"PNG {cname} chunk (length runs past the file)")
                break
            if ctype in PNG_PIXEL_CHUNKS:
                bound = PNG_PIXEL_CHUNKS[ctype]
                if bound and not bound[0] <= n <= bound[1]:
                    found.append(f"PNG {cname} chunk (unexpected length)")
            elif ctype == b"eXIf":
                vet(i, i + 12 + n, exif_fields(body))
            elif ctype == b"iCCP":
                pname, _, rest = body.partition(b"\x00")
                try:
                    prof = zlib.decompress(rest[1:])
                except zlib.error:
                    prof = b""
                fields = icc_fields(prof)
                if not _safe_text(pname.decode("latin-1")):
                    fields.append("PNG iCCP profile name")
                vet(i, i + 12 + n, fields)
            elif ctype == b"iTXt" and body.startswith(b"XML:com.adobe.xmp\x00"):
                rest = body[len(b"XML:com.adobe.xmp\x00"):]
                comp = rest[:1] == b"\x01"
                parts = rest[2:].split(b"\x00", 2)           # language tag, translated keyword, text
                packet = parts[2] if len(parts) == 3 else b""
                if comp:
                    try:
                        packet = zlib.decompress(packet)
                    except zlib.error:
                        packet = b""
                vet(i, i + 12 + n, xmp_fields(packet), XMP_SIGS)
            elif ctype in (b"tEXt", b"zTXt", b"iTXt"):
                found.append(_png_text_keyword(cname, body))
            else:
                found.append(f"PNG {cname} chunk")          # tIME, private chunks, anything unknown
            i += 12 + n
            if ctype == b"IEND":
                if data[i:].strip(b"\x00"):
                    found.append("data after PNG IEND")
                break
    elif kind == "jpeg":
        end = data.rfind(b"\xff\xd9")
        if end != -1 and data[end + 2:].strip(b"\x00\r\n"):
            found.append("data after JPEG end")
        icc_parts = []
        i = 2
        while i < len(data):
            if data[i] != 0xFF:
                # entropy-coded scan data: skip to the next real marker (not stuffing, not RSTn)
                j = data.find(b"\xff", i)
                while j != -1 and j + 1 < len(data) and (data[j + 1] == 0 or 0xD0 <= data[j + 1] <= 0xD7
                                                          or data[j + 1] == 0xFF):
                    j = data.find(b"\xff", j + 2 if data[j + 1] != 0xFF else j + 1)
                if j == -1:
                    break
                i = j
                continue
            while i < len(data) and data[i] == 0xFF:   # fill bytes before a marker
                i += 1
            if i >= len(data):
                break
            marker = data[i]
            i += 1
            if marker == 0xD9:                   # end of image: nothing may follow the FIRST one
                if data[i:].strip(b"\x00\r\n"):
                    found.append("data after JPEG end")
                break
            if 0xD0 <= marker <= 0xD8 or marker == 0x01:
                continue                         # standalone markers, no length
            n = int.from_bytes(data[i:i + 2], "big")
            seg = data[i + 2:i + n]
            span = (i - 2, i + n)
            if n < 2 or i + n > len(data):
                found.append(f"JPEG segment 0x{marker:02X} (length runs past the file)")
                break
            if 0x02 <= marker <= 0xBF or 0xF0 <= marker <= 0xFD:
                found.append(f"JPEG marker 0x{marker:02X} (unknown segment)")
            if marker == 0xE0:
                if not seg.startswith((b"JFIF\x00", b"JFXX\x00")):
                    found.append("JPEG APP0 (not JFIF)")
            elif marker == 0xE1 and seg.startswith(b"Exif\x00\x00"):
                vet(*span, exif_fields(seg[6:]), [b"Exif\x00\x00"])
            elif marker == 0xE1 and seg.startswith(b"http://ns.adobe.com/xap/1.0/\x00"):
                vet(*span, xmp_fields(seg[29:]), XMP_SIGS)
            elif marker == 0xE1:
                found.append("JPEG APP1 (extended XMP / other)")
            elif marker == 0xE2 and seg.startswith(b"ICC_PROFILE\x00"):
                icc_parts.append((seg[12:13], seg[14:]))
            elif marker == 0xE2:
                found.append("JPEG APP2 (MPF / FlashPix / other)")
            elif marker == 0xED and seg.startswith(b"Photoshop 3.0\x00"):
                found.extend(iptc_fields(seg[14:]))
            elif marker == 0xEE:
                if not seg.startswith(b"Adobe"):
                    found.append("JPEG APP14 (not Adobe colour)")
            elif 0xE3 <= marker <= 0xEF:
                found.append(f"JPEG APP{marker - 0xE0} segment")
            elif marker == 0xFE:
                found.append("JPEG comment")
            i += n
        if icc_parts:
            found.extend(icc_fields(b"".join(p for _, p in sorted(icc_parts))))
    elif kind == "webp":
        i = 12
        while i + 8 <= len(data):
            ctype = data[i:i + 4]
            n = int.from_bytes(data[i + 4:i + 8], "little")
            body = data[i + 8:i + 8 + n]
            span = (i, i + 8 + n)
            if ctype in (b"VP8 ", b"VP8L", b"VP8X", b"ALPH", b"ANIM", b"ANMF"):
                pass
            elif ctype == b"EXIF":
                vet(*span, exif_fields(body[6:] if body.startswith(b"Exif\x00\x00") else body),
                    [b"Exif\x00\x00"] if body.startswith(b"Exif\x00\x00") else [])
            elif ctype == b"XMP ":
                vet(*span, xmp_fields(body), XMP_SIGS)
            elif ctype == b"ICCP":
                vet(*span, icc_fields(body))
            else:
                found.append(f"WebP {ctype.decode('latin-1').strip() if ctype.isascii() else ctype.hex()} chunk")
            i += 8 + n + (n & 1)
    elif kind == "gif":
        i = 13
        if len(data) > 10 and data[10] & 0x80:
            i += 3 * (2 ** ((data[10] & 0x07) + 1))
        while i < len(data):
            b = data[i]
            if b == 0x3B:
                if data[i + 1:].strip(b"\x00\r\n"):
                    found.append("data after GIF trailer")
                break
            if b == 0x21 and i + 1 < len(data):
                label = data[i + 1]
                j = i + 2
                first = data[j + 1:j + 1 + data[j]] if j < len(data) else b""
                if label == 0xFE:
                    found.append("GIF comment")
                elif label == 0xFF and not first.startswith((b"NETSCAPE2.0", b"ANIMEXTS1.0")):
                    found.append("GIF application extension (XMP/other)")
                elif label not in (0xF9, 0x01, 0xFF):
                    found.append(f"GIF extension 0x{label:02X}")
                while j < len(data) and data[j]:
                    j += 1 + data[j]
                i = j + 1
            elif b == 0x2C:
                j = i + 10
                if i + 9 < len(data) and data[i + 9] & 0x80:
                    j += 3 * (2 ** ((data[i + 9] & 0x07) + 1))
                j += 1                           # LZW minimum code size
                while j < len(data) and data[j]:
                    j += 1 + data[j]
                i = j + 1
            else:
                found.append("GIF block not understood")
                break
        else:
            found.append("GIF without a trailer")
    if kind:
        budget = [(a, b, list(sigs)) for a, b, sigs in vetted]
        for sig, what in METADATA_SIGNATURES:
            start = data.find(sig)
            while start != -1:
                span = next((v for v in budget if v[0] <= start < v[1] and sig in v[2]), None)
                if span is None:
                    found.append(what + " outside a parsed block")
                    break
                span[2].remove(sig)
                start = data.find(sig, start + 1)
    return list(dict.fromkeys(found))


# ------------------------------------------------------------------ git plumbing
BLOB_SCAN_LIMIT = 25 * 1024 * 1024   # raw-scan binaries up to this size (videos are skipped)


def git(*args, input_bytes=None):
    """Run git; return stdout decoded losslessly, or raise CalledProcessError."""
    r = subprocess.run(["git", *args], input=input_bytes, capture_output=True, check=True)
    return r.stdout.decode("utf-8", "surrogateescape")


def git_bytes(*args, input_bytes=None):
    return subprocess.run(["git", *args], input=input_bytes, capture_output=True, check=True).stdout


def iter_blobs(shas):
    """Yield (sha, bytes) one at a time from a single `git cat-file --batch`, so a push full of
    media is never held in memory at once."""
    if not shas:
        return
    p = subprocess.Popen(["git", "cat-file", "--batch"], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    try:
        for sha in shas:
            p.stdin.write((sha + "\n").encode())
            p.stdin.flush()
            head = p.stdout.readline().split()
            if len(head) < 3:
                continue
            data = p.stdout.read(int(head[2]))
            p.stdout.read(1)                     # trailing newline
            yield sha, data
    finally:
        p.stdin.close()
        p.wait()


# Raw-scanning a binary costs ~200 ms/MB with every pattern; scanning only its printable runs
# (what `strings` shows) costs ~20 ms/MB and finds the same ASCII hits.
_PRINTABLE_RUN = re.compile(rb"[\x20-\x7e]{6,}")
_UTF16LE_RUN = re.compile(rb"(?:[\x20-\x7e]\x00){6,}")
_UTF16BE_RUN = re.compile(rb"(?:\x00[\x20-\x7e]){6,}")


def binary_text(data: bytes) -> str:
    runs = [r.decode("ascii") for r in _PRINTABLE_RUN.findall(data)]
    runs += [r.decode("utf-16-le") for r in _UTF16LE_RUN.findall(data)]
    runs += [r.decode("utf-16-be") for r in _UTF16BE_RUN.findall(data)]
    return "\n".join(runs)


def ignored_but_tracked(paths):
    """Subset of `paths` a .gitignore rule matches although git still tracks them.

    gitignore never untracks a file added before the rule, so a hit can sit in a file
    everyone believes is ignored (2026-09-29: 21 hits in 118 MB of ignored-but-tracked
    scaffolding, found only by hand). Best effort — returns nothing outside a git repo.
    """
    paths = sorted(set(p for p in paths if p))
    if not paths:
        return set()
    try:
        out = subprocess.run(["git", "check-ignore", "--no-index", "--stdin", "-z"],
                             input="\0".join(paths) + "\0", capture_output=True, text=True)
    except OSError:
        return set()
    return set(p for p in out.stdout.split("\0") if p)


ZERO = "0" * 40


def scan_push(local, remote, remote_name, identity):
    """Every check for one pushed ref. Returns a list of (section, line) findings, or raises."""
    if remote == ZERO:
        # new ref: everything not already on THIS remote (another remote may be private)
        rng = [local, "--not", f"--remotes={remote_name}"]
    else:
        rng = [f"{remote}..{local}"]
    git("rev-list", *rng)                                    # fail closed if uncomputable
    findings = []

    # 1. content of every commit in the range
    # every output-shaping option pinned, so no user config (color.ui, diff.noprefix,
    # log.showRoot, diff.external, ...) can blind or misdirect the parser
    patch = git("-c", "color.ui=never", "log", "-p", "--root", "--no-color", "--no-ext-diff",
                "--no-textconv", "--src-prefix=a/", "--dst-prefix=b/", "--no-renames", "--unified=0",
                "--diff-merges=first-parent", "--format=%x00commit %h", *rng)
    added, binary_paths = parse_patch(patch)
    content = []
    seen = set()
    for commit, path, lineno, line in added:
        for label, sample in match_line(line, identity):
            key = (path, line.strip(), label)
            if key not in seen:
                seen.add(key)
                content.append((path, lineno, commit, label, sample))
    stale = ignored_but_tracked(c[0] for c in content)
    for path, lineno, commit, label, sample in content:
        mark = "  [gitignored but still tracked]" if path in stale else ""
        findings.append(("content", f"{path}:{lineno} @{commit}: {label}: {sample}{mark}"))
    if stale:
        findings.append(("stale", ""))

    # 2. new objects: file names, binary content, image metadata
    objects = []
    for raw in git("rev-list", "--objects", *rng).split("\n"):
        sha, _, path = raw.partition(" ")
        if path:
            objects.append((sha, path))
    names = git("log", "--root", "--no-renames", "--diff-merges=first-parent", "--format=",
                "--name-only", "-z", *rng)
    for path in sorted(set(n.strip("\n") for n in names.split("\0") if n.strip("\n"))):
        for label, sample in match_line(path, identity):
            findings.append(("content", f"{path!r}: file NAME matches {label}: {sample}"))
    checks = git_bytes("cat-file", "--batch-check=%(objectname) %(objecttype) %(objectsize)",
                       input_bytes=("\n".join(s for s, _ in objects) + "\n").encode()).decode().split("\n")
    sizes = {}
    for c in checks:
        parts = c.split()
        if len(parts) == 3 and parts[1] == "blob":
            sizes[parts[0]] = int(parts[2])
    wanted = [s for s, _ in objects if s in sizes and sizes[s] <= BLOB_SCAN_LIMIT]
    paths_of = {}
    for sha, path in objects:
        paths_of.setdefault(sha, []).append(path)
    for sha, data in iter_blobs(list(dict.fromkeys(wanted))):
        for path in paths_of[sha]:
            if path in SELF_FILES:
                continue
            if b"\0" in data[:8000] or path in binary_paths:
                for label, sample in set(match_line(binary_text(data), identity, binary=True)):
                    findings.append(("content", f"{path} (binary): {label}: {sample}"))
            if identity:
                for what in image_metadata(data, path):
                    findings.append(("image", f"{path}: {what}"))
    return findings + scan_messages_and_authors(local, rng, identity)


def scan_messages_and_authors(local, rng, identity):
    findings = []

    # 3. commit and tag messages (+ tagger identity)
    messages = git("log", "--format=%x00%h%n%B", *rng)
    obj = local
    for _ in range(20):                          # a tag can point at a tag: peel the chain
        if git("cat-file", "-t", obj).strip() != "tag":
            break
        body = git("cat-file", "tag", obj)
        messages += "\x00tag\n" + body
        obj = body.split("\n", 1)[0].split()[1]
    ref = ""
    for line in messages.split("\n"):
        if line.startswith("\x00"):
            ref = line[1:]
            continue
        for label, sample in match_line(line, identity):
            findings.append(("content", f"{ref} message: {label}: {sample}"))
        if identity and ref == "tag" and line.startswith("tagger "):
            email = line[line.find("<") + 1:line.find(">")]
            if not NOREPLY_RE.search(email):
                findings.append(("author", f"tag tagger: {redact(email) if email else '(empty)'}"))

    # 4. authors
    if identity:
        log = git("log", "--format=%h%x09%ae%x09%ce", *rng)
        for sha, role, email in scan_authors(log):
            findings.append(("author", f"{sha} {role}: {email}"))
        for line in git("log", "--format=%h%x09%an%x09%cn", *rng).split("\n"):
            sha, _, names = line.partition("\t")
            for label, sample in match_line(names, identity):
                findings.append(("author", f"{sha} author/committer NAME: {label}: {sample}"))
    return findings


# ------------------------------------------------------------------ main
BYPASS = "If every hit above is a false positive, bypass with SKIP_SECRET_SCAN=1 git push."


def err(*lines):
    for line in lines:
        print(line, file=sys.stderr)


def report(findings, ref_desc):
    content = [f for s, f in findings if s == "content"]
    authors = [f for s, f in findings if s == "author"]
    images = [f for s, f in findings if s == "image"]
    err("", f"pre-push: BLOCKED — {ref_desc}")
    if content:
        err("The outgoing commits contain what looks like a credential, a private path/host"
            " or a personal identifier:")
        err(*(f"  {f}" for f in content))
    if any(s == "stale" for s, _ in findings):
        err("Some hits are in files a .gitignore rule matches but git still tracks (gitignore",
            "never untracks a file added before the rule). List them all with",
            "    git ls-files -ci --exclude-standard",
            "and untrack with git rm --cached <path> (the files stay on disk).")
    if authors:
        err("Outgoing commits/tags carry a non-noreply email (commit metadata is permanent once pushed):")
        err(*(f"  {f}" for f in authors))
        err("Fix: git config user.email '<id>+<user>@users.noreply.github.com', then rewrite the",
            "unpushed commits (git rebase -r <base> --exec 'git commit --amend --no-edit --reset-author')",
            "and purge the reflog afterwards, per the checklist.")
    if images:
        err("Outgoing images carry metadata that is private or not known to be harmless",
            "(field names only; values are never printed):")
        err(*(f"  {f}" for f in images))
        err("Open and look at each image, then strip it (mogrify -strip <file>, or exiftool -all= <file>)",
            "and amend the commit that added it. Pixel size, resolution, orientation, colour profile and",
            "a 'Screenshot' comment are allowed; a field named here is not.")
    err("", BYPASS)


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    identity = "--identity" in argv

    if "--push" in argv:
        i = argv.index("--push")
        local, remote, remote_name = argv[i + 1:i + 4]
        try:
            findings = scan_push(local, remote, remote_name, identity)
        except (subprocess.CalledProcessError, OSError, ValueError) as e:
            detail = getattr(e, "stderr", b"") or b""
            err("", f"pre-push: BLOCKED — could not compute what {local[:12]} would send "
                f"({detail.decode(errors='replace').strip()[:200] or e}).",
                "Is the remote ref fetched? Try 'git fetch' and push again.")
            return 1
        if findings:
            report(findings, f"pushing {local[:12]}")
            return 1
        return 0

    # stdin mode: a unified diff (kept for ad-hoc use and the unit tests)
    findings = scan(sys.stdin.buffer.read().decode("utf-8", "surrogateescape"), identity=identity)
    if not findings:
        return 0
    report([("content", f"{p}:{n}: {label}: {s}") for p, n, label, s in findings], "the diff")
    return 1


if __name__ == "__main__":
    sys.exit(main())
