#!/usr/bin/env python3
"""
Regression test for githooks/scan_diff.py and this repo's githooks/pre-push.

    python3 githooks/test_scan_diff.py

Byte-identical in every repo that carries the scanner. The scanner unit tests run everywhere;
the end-to-end cases build a throwaway repo + bare remote, copy THIS repo's githooks/ into it
and push for real, so the hook's wiring is tested, not a stand-in. Which end-to-end cases
apply is read from the hook itself: a hook that passes `--identity` (the public_tools repos)
must block the identity grep, non-noreply authors and image metadata; one that does not (the
website, escape_rooms, gcms) must let a contact address through.

FAILURE MODES PINNED
- 2026-09-29: the scan blocked a website push with 21 redacted hits and NO file or line, so
  finding them took a round of manual git archaeology on the Mac — every hit is now `path:line`,
  flagged when the file is gitignored-but-still-tracked (the actual cause that day).
- 2026-09-28: public_tools' pre-publish checklist (identity grep, noreply author, stripped
  images) was a manual step with no gate; the two covers shipped as "stripped" still carried
  eXIf + XMP chunks. `--identity` now enforces the grep, noreply authors and image metadata.
- 2026-10-08 adversarial check of the first version: a secret added then deleted inside the push
  range, a binary or `-diff` file, a commit message, an added line reading `++ githooks/x`
  (parsed as a file header and skipped), a U+2028 in a line, `color.ui=always`,
  `/home/bustalab_old`, `/HOME/BUSTALAB`, a `src/githooks/` path, a missing scanner, a padded
  JPEG, a GIF comment and a TIFF renamed .png all got through. Each has a case below.
- 2026-10-08 second adversarial pass: `log.showRoot=false` (root commit unscanned),
  `diff.noprefix` (a `b/githooks/...` path posing as the scanner), `.env`-style unquoted tokens,
  lowercase `api_key = "..."`, Canvas `NNNN~...` and Bearer tokens, a second assignment on a line
  after a placeholder, `labbust0037`, an x500 author NAME, a tag pointing at a tag, a new file
  whose blob was already on the remote, UTF-16 text, JPEG data after the scan, a private PNG
  chunk, AVIF `avis` and BigTIFF all got through. Each has a case below.
- 2026-10-08 (afternoon): the blanket "no EXIF/XMP anywhere" rule blocked the two published
  covers, whose metadata is macOS screenshot data (pixel size, "Screenshot", a display colour
  profile). Lucas: harmless metadata must not block; private metadata must. Images are now parsed
  against an allowlist — the covers' own metadata is replicated in `screenshot_png()` and must
  pass, and GPS, make/model, artist, a non-Screenshot comment, an XMP edit history, an unknown EXIF
  tag, an EXIF thumbnail, IPTC, a display serial and a person's name in an ICC copyright must each
  block and be NAMED, with the value never printed.
"""
import importlib.util
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
import zlib

sys.dont_write_bytecode = True          # never leave a __pycache__ in a published repo
HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("scan_diff", os.path.join(HERE, "scan_diff.py"))
scan_diff = importlib.util.module_from_spec(spec)
spec.loader.exec_module(scan_diff)

HOOK = os.path.join(HERE, "pre-push")
HOOK_TEXT = open(HOOK, encoding="utf-8").read() if os.path.exists(HOOK) else ""
# read from the hook's CODE, not its comments (the website hook's header mentions the flag)
IDENTITY_HOOK = any("--identity" in l for l in HOOK_TEXT.splitlines() if not l.lstrip().startswith("#"))

NOREPLY = "1+test@users.noreply.github.com"
ENV = dict(os.environ, GIT_AUTHOR_NAME="t", GIT_COMMITTER_NAME="t",
           GIT_AUTHOR_EMAIL=NOREPLY, GIT_COMMITTER_EMAIL=NOREPLY,
           GIT_CONFIG_GLOBAL=os.devnull, GIT_CONFIG_NOSYSTEM="1")
ENV.pop("SKIP_SECRET_SCAN", None)

# fragments assembled at runtime so this file never contains a literal the scanner flags
HOME = "/home/" + "bustalab"
UMN_MAIL = "someone" + "@d.umn" + ".edu"
UMN_LINK = "https://chem.d.umn" + ".edu/people"
GMAIL = "someone" + "@gmail" + ".com"
TOKEN = "sk-" + "FAKE1234567890abcdefFAKE"


def png(chunks=()):
    def chunk(t, d):
        return len(d).to_bytes(4, "big") + t + d + zlib.crc32(t + d).to_bytes(4, "big")
    ihdr = (1).to_bytes(4, "big") * 2 + bytes([8, 0, 0, 0, 0])
    body = chunk(b"IHDR", ihdr) + b"".join(chunk(t, d) for t, d in chunks)
    return b"\x89PNG\r\n\x1a\n" + body + chunk(b"IDAT", zlib.compress(b"\x00\x00")) + chunk(b"IEND", b"")


def tiff(ifd0, exif=None, gps=False, thumb=False):
    """A big-endian TIFF/EXIF block. ifd0/exif: list of (tag, type, count, value-bytes)."""
    def ifd(entries, base, nxt=0):
        entries = sorted(entries)
        head = len(entries).to_bytes(2, "big")
        extra, body = b"", b""
        data_at = base + 2 + 12 * len(entries) + 4
        for tag, typ, count, val in entries:
            if len(val) <= 4:
                field = val.ljust(4, b"\x00")
            else:
                field = (data_at + len(extra)).to_bytes(4, "big")
                extra += val + (b"\x00" if len(val) & 1 else b"")
            body += tag.to_bytes(2, "big") + typ.to_bytes(2, "big") + count.to_bytes(4, "big") + field
        return head + body + nxt.to_bytes(4, "big") + extra

    def size(entries):
        return len(ifd(entries, 0))

    ifd0 = list(ifd0)
    ptrs = []
    if exif is not None:
        ptrs.append(0x8769)
    if gps:
        ptrs.append(0x8825)
    for t in ptrs:
        ifd0.append((t, 4, 1, b"\x00\x00\x00\x00"))
    base0 = 8
    after0 = base0 + size(ifd0)
    sub = {}
    cursor = after0
    if exif is not None:
        sub[0x8769] = (cursor, exif)
        cursor += size(exif)
    if gps:
        gps_entries = [(0x0002, 5, 3, b"\x00\x00\x00\x2e\x00\x00\x00\x01" * 3)]   # latitude
        sub[0x8825] = (cursor, gps_entries)
        cursor += size(gps_entries)
    thumb_at = cursor if thumb else 0
    ifd0 = [(t, ty, c, sub[t][0].to_bytes(4, "big") if t in sub else v) for t, ty, c, v in ifd0]
    out = b"MM\x00*" + base0.to_bytes(4, "big") + ifd(ifd0, base0, thumb_at)
    for t in sorted(sub):
        out += ifd(sub[t][1], sub[t][0])
    if thumb:
        out += ifd([(0x0103, 3, 1, b"\x00\x06")], thumb_at)
    return out


def short(v):
    return v.to_bytes(2, "big")


def long_(v):
    return v.to_bytes(4, "big")


RES = (1).to_bytes(4, "big") * 0 + (144).to_bytes(4, "big") + (1).to_bytes(4, "big")
SCREENSHOT_COMMENT = b"ASCII\x00\x00\x00Screenshot"
SCREENSHOT_IFD0 = [(0x011A, 5, 1, RES), (0x011B, 5, 1, RES), (0x0128, 3, 1, short(2))]
SCREENSHOT_EXIF = [(0x9286, 7, len(SCREENSHOT_COMMENT), SCREENSHOT_COMMENT),
                   (0xA002, 4, 1, long_(2038)), (0xA003, 4, 1, long_(598))]


def xmp(props):
    return ('<x:xmpmeta xmlns:x="adobe:ns:meta/" x:xmptk="XMP Core 6.0.0"><rdf:RDF '
            'xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#"><rdf:Description rdf:about="" '
            'xmlns:exif="http://ns.adobe.com/exif/1.0/" xmlns:xmpMM="http://ns.adobe.com/xap/1.0/mm/" '
            'xmlns:stEvt="http://ns.adobe.com/xap/1.0/sType/ResourceEvent#" '
            'xmlns:xmp="http://ns.adobe.com/xap/1.0/">' + props +
            '</rdf:Description></rdf:RDF></x:xmpmeta>').encode()


SCREENSHOT_XMP = ("<exif:PixelYDimension>598</exif:PixelYDimension>"
                  "<exif:PixelXDimension>2038</exif:PixelXDimension>"
                  "<exif:UserComment>Screenshot</exif:UserComment>")


def icc(desc="Display", cprt="Copyright Apple Inc., 2026", serial=0xFD626D62, extra=(), desc_raw=None):
    """A minimal display profile shaped like macOS's: desc, cprt, wtpt, mmod (+ extra tags)."""
    def text_tag(t):
        return b"text\x00\x00\x00\x00" + t.encode() + b"\x00"
    def desc_tag(t):
        return b"desc\x00\x00\x00\x00" + (len(t) + 1).to_bytes(4, "big") + t.encode() + b"\x00" * 80
    mmod = b"mmod" + b"\x00" * 4 + long_(0x0610) + long_(0xA061) + long_(serial) + long_(0) + b"\x00" * 16
    tags = [(b"desc", desc_raw or desc_tag(desc)), (b"cprt", text_tag(cprt)), (b"wtpt", b"XYZ " + b"\x00" * 16),
            (b"mmod", mmod)] + list(extra)
    table, blob = b"", b""
    off = 132 + 12 * len(tags)
    for sig, data in tags:
        data += b"\x00" * (-len(data) % 4)
        table += sig + long_(off + len(blob)) + long_(len(data))
        blob += data
    header = bytearray(128)
    header[4:8], header[12:20], header[36:40] = b"appl", b"mntrRGB ", b"acsp"
    prof = bytes(header) + long_(len(tags)) + table + blob
    return long_(len(prof)) + prof[4:]


def iccp(profile):
    return (b"iCCP", b"ICC Profile\x00\x00" + zlib.compress(profile))


def screenshot_png(**over):
    """A PNG carrying exactly the metadata the two published covers carry (2026-10-08)."""
    parts = {"iccp": iccp(icc()), "exif": (b"eXIf", tiff(SCREENSHOT_IFD0, SCREENSHOT_EXIF)),
             "phys": (b"pHYs", b"\x00\x00\x16%\x00\x00\x16%\x01"),
             "xmp": (b"iTXt", b"XML:com.adobe.xmp\x00\x00\x00\x00\x00" + xmp(SCREENSHOT_XMP)),
             "idot": (b"iDOT", b"\x00" * 28)}
    parts.update(over)
    return png([c for c in parts.values() if c])


def find_covers():
    """The real published covers, wherever this copy of the test sits relative to them."""
    out = []
    for rel in ("..", os.path.join("..", "..", "doi2abstracts"), os.path.join("..", "..", "canvas-assignments"),
                os.path.join("..", "..", "..", "doi2abstracts"), os.path.join("..", "..", "..", "canvas-assignments")):
        p = os.path.normpath(os.path.join(HERE, rel, "cover.png"))
        if os.path.exists(p) and p not in out:
            out.append(p)
    return out


class ImageAllowlist(unittest.TestCase):
    def assertBlocks(self, data, field, secret=None):
        found = scan_diff.image_metadata(data, "x")
        self.assertTrue(any(field in f for f in found), f"{field!r} not named in {found}")
        if secret:
            self.assertNotIn(secret, " ".join(found), "a private VALUE was printed")

    def test_screenshot_metadata_passes(self):
        self.assertEqual(scan_diff.image_metadata(screenshot_png(), "cover.png"), [])

    @unittest.skipUnless(find_covers(), "the published covers are not beside this copy")
    def test_real_published_covers_pass(self):
        for p in find_covers():
            with open(p, "rb") as f:
                self.assertEqual(scan_diff.image_metadata(f.read(), p), [], p)

    def test_gps_blocks(self):
        self.assertBlocks(png([(b"eXIf", tiff(SCREENSHOT_IFD0, SCREENSHOT_EXIF, gps=True))]), "GPSInfo")

    def test_make_and_model_block(self):
        ifd0 = SCREENSHOT_IFD0 + [(0x010F, 2, 6, b"Apple\x00"), (0x0110, 2, 10, b"iPhone 15\x00")]
        found = scan_diff.image_metadata(png([(b"eXIf", tiff(ifd0, SCREENSHOT_EXIF))]), "x")
        self.assertIn("EXIF Make", found)
        self.assertIn("EXIF Model", found)
        self.assertNotIn("iPhone", " ".join(found))

    def test_artist_blocks(self):
        ifd0 = SCREENSHOT_IFD0 + [(0x013B, 2, 12, b"Jane Doe X1\x00")]
        self.assertBlocks(png([(b"eXIf", tiff(ifd0, SCREENSHOT_EXIF))]), "EXIF Artist", "Jane")

    def test_non_screenshot_user_comment_blocks(self):
        c = b"ASCII\x00\x00\x00Shot at home"
        exif = [(0x9286, 7, len(c), c)]
        self.assertBlocks(png([(b"eXIf", tiff(SCREENSHOT_IFD0, exif))]), "EXIF UserComment", "home")

    def test_unknown_exif_tag_blocks(self):
        ifd0 = SCREENSHOT_IFD0 + [(0x7777, 3, 1, short(1))]
        self.assertBlocks(png([(b"eXIf", tiff(ifd0, SCREENSHOT_EXIF))]), "EXIF tag 0x7777")

    def test_serial_owner_makernote_thumbnail_block(self):
        exif = SCREENSHOT_EXIF + [(0xA431, 2, 8, b"SN12345\x00"), (0xA430, 2, 4, b"Bob\x00"),
                                  (0x927C, 7, 4, b"abcd")]
        found = scan_diff.image_metadata(png([(b"eXIf", tiff(SCREENSHOT_IFD0, exif, thumb=True))]), "x")
        for name in ("BodySerialNumber", "CameraOwnerName", "MakerNote", "IFD1 (thumbnail"):
            self.assertTrue(any(name in f for f in found), (name, found))

    def test_software_with_a_path_blocks_plain_software_passes(self):
        ok = SCREENSHOT_IFD0 + [(0x0131, 2, 10, b"GIMP 2.10\x00")]
        self.assertEqual(scan_diff.image_metadata(png([(b"eXIf", tiff(ok, SCREENSHOT_EXIF))]), "x"), [])
        bad = SCREENSHOT_IFD0 + [(0x0131, 2, 14, b"/opt/jane/bin\x00")]
        self.assertBlocks(png([(b"eXIf", tiff(bad, SCREENSHOT_EXIF))]), "EXIF Software", "jane")

    def test_xmp_history_trail_blocks(self):
        trail = ('<xmpMM:History><rdf:Seq><rdf:li rdf:parseType="Resource"><stEvt:action>saved'
                 '</stEvt:action><stEvt:softwareAgent>Jane Laptop</stEvt:softwareAgent></rdf:li>'
                 '</rdf:Seq></xmpMM:History><xmpMM:DocumentID>xmp.did:1234</xmpMM:DocumentID>')
        data = screenshot_png(xmp=(b"iTXt", b"XML:com.adobe.xmp\x00\x00\x00\x00\x00" + xmp(SCREENSHOT_XMP + trail)))
        found = scan_diff.image_metadata(data, "x")
        self.assertIn("XMP xmpMM:History", found)
        self.assertIn("XMP xmpMM:DocumentID", found)
        self.assertNotIn("Jane", " ".join(found))

    def test_xmp_creator_tool_with_path_and_gps_block(self):
        props = '<xmp:CreatorTool>C:\\Users\\jane\\tool.exe</xmp:CreatorTool><exif:GPSLatitude>46,47N</exif:GPSLatitude>'
        found = scan_diff.image_metadata(png([(b"iTXt", b"XML:com.adobe.xmp\x00\x00\x00\x00\x00" + xmp(props))]), "x")
        self.assertIn("XMP xmp:CreatorTool", found)
        self.assertIn("XMP exif:GPSLatitude", found)

    def test_xmp_non_screenshot_comment_blocks(self):
        data = png([(b"iTXt", b"XML:com.adobe.xmp\x00\x00\x00\x00\x00" + xmp("<exif:UserComment>Lab 12</exif:UserComment>"))])
        self.assertBlocks(data, "XMP exif:UserComment", "Lab 12")

    def test_iptc_blocks_and_names_the_dataset(self):
        iptc = b"\x1c\x02\x50\x00\x08Jane Doe" + b"\x1c\x02\x5a\x00\x06Duluth"
        irb = b"8BIM\x04\x04\x00\x00" + long_(len(iptc)) + iptc
        seg = b"Photoshop 3.0\x00" + irb
        jpeg = b"\xff\xd8\xff\xed" + short(len(seg) + 2) + seg + b"\xff\xd9"
        found = scan_diff.image_metadata(jpeg, "x")
        self.assertIn("IPTC 2:80 By-line", found)
        self.assertIn("IPTC 2:90 City", found)
        self.assertNotIn("Duluth", " ".join(found))

    def test_icc_display_serial_and_person_copyright_block(self):
        self.assertBlocks(screenshot_png(iccp=iccp(icc(serial=0x12345678))), "ICC mmod")
        self.assertBlocks(screenshot_png(iccp=iccp(icc(cprt="Copyright Jane Doe 2026"))), "ICC cprt", "Jane")
        self.assertBlocks(screenshot_png(iccp=iccp(icc(desc="/Users/jane/x.icc"))), "ICC desc", "jane")
        self.assertBlocks(screenshot_png(iccp=iccp(icc(extra=[(b"dmdd", b"x" * 12)]))), "ICC dmdd")
        self.assertBlocks(screenshot_png(iccp=iccp(icc(extra=[(b"zzzz", b"x" * 12)]))), "ICC zzzz")

    def test_software_naming_a_person_or_home_blocks(self):
        for sw in ("Jane's MacBook Pro", "home tool 1.0", "Users tool"):
            ifd0 = SCREENSHOT_IFD0 + [(0x0131, 2, len(sw) + 1, sw.encode() + b"\x00")]
            self.assertBlocks(png([(b"eXIf", tiff(ifd0, SCREENSHOT_EXIF))]), "EXIF Software", sw)

    def test_pixel_chunk_with_room_for_text_blocks(self):
        self.assertBlocks(png([(b"pHYs", b"\x00" * 9 + b"lat 46.78N lon 92.1W")]), "PNG pHYs chunk (unexpected length)")
        self.assertBlocks(png([(b"iDOT", b"\x00" * 28 + b"Jane")]), "PNG iDOT chunk (unexpected length)")

    def test_non_iptc_photoshop_resource_is_named(self):
        irb = b"8BIM\x04\x0c\x00\x00" + long_(4) + b"abcd"           # 0x040C: thumbnail
        seg = b"Photoshop 3.0\x00" + irb
        jpeg = b"\xff\xd8\xff\xed" + short(len(seg) + 2) + seg + b"\xff\xd9"
        self.assertBlocks(jpeg, "Photoshop resource 0x040C")

    def test_adversarial_round_one(self):
        """Each of these got past the first allowlist (Sonnet adversarial check, 2026-10-08)."""
        def sw(text):
            ifd0 = SCREENSHOT_IFD0 + [(0x0131, 2, len(text), text)]
            return png([(b"eXIf", tiff(ifd0, SCREENSHOT_EXIF))])
        self.assertBlocks(sw(b"Photoshop\x00 Jane Doe\x00"), "EXIF Software", "Jane")    # after the NUL
        self.assertBlocks(sw(b"janes-mbp.local\x00"), "EXIF Software", "janes")         # a hostname
        self.assertBlocks(sw(b"Edited by Jane Doe\x00"), "EXIF Software", "Jane")       # an attribution
        self.assertBlocks(screenshot_png(iccp=iccp(icc(cprt="Copyright Apple Inc.\x00Jane Doe"))), "ICC cprt")
        uni = "C:/Users/jane/cal.icc".encode("utf-16-be")                 # the UTF-16 half of `desc`
        desc = b"desc\x00\x00\x00\x00" + long_(8) + b"Display\x00" + long_(0) + long_(len(uni) // 2) + uni + b"\x00" * 70
        self.assertBlocks(screenshot_png(iccp=iccp(icc(desc_raw=desc))), "ICC desc", "Jane")
        xp = lambda body: png([(b"iTXt", b"XML:com.adobe.xmp\x00\x00\x00\x00\x00" + body)])
        ok = xmp("<exif:PixelXDimension>10</exif:PixelXDimension>")
        self.assertEqual(scan_diff.image_metadata(xp(ok), "x"), [])
        self.assertBlocks(xp(ok + xmp("<xmpMM:DocumentID>x</xmpMM:DocumentID>")), "XMP packet")
        self.assertBlocks(xp(ok.replace(b"<exif:PixelX", b"<!-- Jane --><exif:PixelX")), "XMP comment", "Jane")
        self.assertBlocks(xp(ok.replace(b"<exif:PixelXDimension>", b'<exif:PixelXDimension xmp:Author="Jane">')),
                          "XMP exif:PixelXDimension (attributes or stray text)", "Jane")
        self.assertBlocks(xp(ok.replace(b"</rdf:Description>", b"</rdf:Description>Jane Doe")), "XMP", "Jane")
        self.assertBlocks(xp(ok.replace(b"XMP Core 6.0.0", b"Jane Doe laptop")), "XMP xmptk attribute", "Jane")
        second_rdf = ok.replace(b"</x:xmpmeta>", b'<rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#"/></x:xmpmeta>')
        self.assertBlocks(xp(second_rdf), "XMP packet (not exactly one rdf:RDF)")
        self.assertBlocks(b"8BPS" + b"\x00" * 30, "cannot be checked")                      # a PSD
        self.assertBlocks(b"\x00\x00\x00\x18ftypcrx " + b"\x00" * 12, "cannot be checked")  # Canon CR3
        clean_gif = b"GIF89a" + b"\x01\x00\x01\x00\x00\x00\x00" + b"\x3b"
        self.assertBlocks(clean_gif + b"GPS Jane", "data after GIF trailer", "Jane")
        self.assertBlocks(b"\xff\xd8\xff\xdb\xff\xf0" + b"\x00" * 10 + b"\xff\xfe\x00\x06Jane\xff\xd9",
                          "length runs past the file")
        self.assertBlocks(b"\xff\xd8\xff\xf5\x00\x06Jane\xff\xd9", "JPEG marker 0xF5", "Jane")
        hidden = png()[:-12] + (10 ** 6).to_bytes(4, "big") + b"IDAT"
        self.assertBlocks(hidden + b"\x00\x00\x00\x0atEXtAuthor\x00Jane", "length runs past the file", "Jane")
        self.assertBlocks(screenshot_png(iccp=iccp(icc(cprt="Copyright Apple Inc. 12 Elm St"))), "ICC cprt")
        self.assertBlocks(b"\xff\xd8\xff\xd9Jane Doe\xff\xd9", "data after JPEG end", "Jane")
        # a second EXIF block hidden in the slack of an allowed one: its parse accounts for one only
        seg = b"Exif\x00\x00" + tiff(SCREENSHOT_IFD0, SCREENSHOT_EXIF) + b"Exif\x00\x00" + tiff([(0x8825, 4, 1, long_(8))])
        self.assertBlocks(b"\xff\xd8\xff\xe1" + short(len(seg) + 2) + seg + b"\xff\xd9", "EXIF block outside a parsed block")
        # a second XMP signature inside a vetted span is not covered by that span's parse
        jicc_seg = b"ICC_PROFILE\x00\x01\x01" + icc(extra=[(b"vcgt", b"vcgt\x00\x00\x00\x00<x:xmpmeta")])
        jicc = b"\xff\xd8\xff\xe2" + short(len(jicc_seg) + 2) + jicc_seg + b"\xff\xd9"
        self.assertBlocks(jicc, "XMP packet outside a parsed block")

    def test_adversarial_round_two(self):
        xp = lambda body: png([(b"iTXt", b"XML:com.adobe.xmp\x00\x00\x00\x00\x00" + body)])
        ok = xmp("<exif:PixelXDimension>10</exif:PixelXDimension>")
        self.assertBlocks(xp(ok.replace(b"XMP Core 6.0.0", b"XMP Core jane doe 12 elm st")), "XMP xmptk attribute", "jane")
        self.assertBlocks(xp(ok.replace(b"XMP Core 6.0.0", b"XMP Core 6.0.0 Jane Doe")), "XMP xmptk attribute", "Jane")
        self.assertBlocks(xp(ok.replace(b"<exif:PixelX", b"<?foo Jane?><exif:PixelX")), "XMP processing instruction", "Jane")
        self.assertBlocks(xp(ok.replace(b"<exif:PixelXDimension>", b'<exif:PixelXDimension xml:lang="jane@home">')),
                          "XMP exif:PixelXDimension", "jane")
        self.assertBlocks(xp(ok.replace(b'rdf:about=""', b'rdf:about="" xmlns:q="http://host.local/jane"')),
                          "XMP (unrecognised namespace declared)", "jane")
        self.assertBlocks(b"\x00\x00\x00\x18ftypiso8\x00\x00\x00\x00heicmif1", "cannot be checked")
        gif = b"GIF89a" + b"\x01\x00\x01\x00\x00\x00\x00"
        self.assertBlocks(gif + b"Jane Doe 12 Elm St", "GIF block not understood", "Jane")
        # a lone unclosed-tag flood must not hang the push (the old lazy regex was quadratic)
        import time
        t0 = time.time()
        scan_diff.xmp_fields(b"<x:xmpmeta " * 40000)
        self.assertLess(time.time() - t0, 2)

    def test_real_world_strings_are_not_false_blocks(self):
        for sw in (b"Pixelmator Pro 3.5 for Mac", b"Paint.NET v5.0", b"Adobe Photoshop 25.0 (Macintosh)", b"GIMP 2.10.36"):
            ifd0 = SCREENSHOT_IFD0 + [(0x0131, 2, len(sw) + 1, sw + b"\x00")]
            self.assertEqual(scan_diff.image_metadata(png([(b"eXIf", tiff(ifd0, SCREENSHOT_EXIF))]), "x"), [], sw)
        for d, c in (("sRGB IEC61966-2.1", "Copyright (c) 1998 Hewlett-Packard Company"),
                     ("Adobe RGB (1998)", "Copyright 1999 Adobe Systems Incorporated"),
                     ("sRGB IEC61966-2.1 : v2", "Copyright 2007 Apple Inc., all rights reserved."),
                     ("Google/Skia/" + "0" * 32, "Google Inc. 2016"), ("Display P3", "Copyright (c) 2015 Canon Inc.")):
            self.assertEqual(scan_diff.image_metadata(screenshot_png(iccp=iccp(icc(desc=d, cprt=c))), "x"), [], (d, c))
        xmptk = b"Adobe XMP Core 9.1-c002 79.dba3da3b5, 2023/12/15-10:42:37        "
        ok = xmp("<exif:PixelXDimension>10</exif:PixelXDimension>").replace(b"XMP Core 6.0.0", xmptk)
        wrapped = b'<?xpacket begin="\xef\xbb\xbf" id="W5M0MpCehiHzreSzNTczkc9d"?>' + ok + b'<?xpacket end="w"?>'
        self.assertEqual(scan_diff.image_metadata(png([(b"iTXt", b"XML:com.adobe.xmp\x00\x00\x00\x00\x00" + wrapped)]), "x"), [])

    def test_png_text_and_time_chunks_block(self):
        self.assertBlocks(png([(b"tEXt", b"Author\x00Jane Doe")]), "PNG tEXt 'Author'", "Jane")
        self.assertBlocks(png([(b"tIME", b"\x07\xea\x0a\x08\x0c\x00\x00")]), "PNG tIME")

    def test_exif_hidden_outside_a_parsed_block(self):
        hidden = b"Exif\x00\x00" + tiff(SCREENSHOT_IFD0 + [(0x013B, 2, 4, b"Bob\x00")], SCREENSHOT_EXIF)
        self.assertBlocks(png() + hidden, "data after PNG IEND")
        self.assertBlocks(png([(b"tRNS", hidden[:200])]), "EXIF block outside a parsed block")

    def test_jpeg_exif_allowlist(self):
        seg = b"Exif\x00\x00" + tiff(SCREENSHOT_IFD0, SCREENSHOT_EXIF)
        ok = b"\xff\xd8\xff\xe1" + short(len(seg) + 2) + seg + b"\xff\xda\x00\x02\x11\x22\xff\xd9"
        self.assertEqual(scan_diff.image_metadata(ok, "x"), [])
        seg = b"Exif\x00\x00" + tiff(SCREENSHOT_IFD0, SCREENSHOT_EXIF, gps=True)
        bad = b"\xff\xd8\xff\xe1" + short(len(seg) + 2) + seg + b"\xff\xd9"
        self.assertBlocks(bad, "GPSInfo")
        # an APP1 after the first scan (progressive JPEG) is still walked
        late = b"\xff\xd8\xff\xda\x00\x02\x11\x22\xff\xe1" + short(len(seg) + 2) + seg + b"\xff\xd9"
        self.assertBlocks(late, "GPSInfo")


class ScannerUnit(unittest.TestCase):
    def test_hit_reports_path_and_line(self):
        diff = ("diff --git a/src/x.py b/src/x.py\n--- a/src/x.py\n+++ b/src/x.py\n"
                "@@ -10,0 +11,3 @@\n+ok = 1\n+ok = 2\n+p = '" + HOME + "/x'\n")
        hits = scan_diff.scan(diff)
        self.assertEqual([(h[0], h[1]) for h in hits], [("src/x.py", 13)])

    def test_line_numbers_track_across_hunks(self):
        diff = ("diff --git a/a.txt b/a.txt\n--- a/a.txt\n+++ b/a.txt\n@@ -1 +1 @@\n-old\n+" + HOME + "\n"
                "@@ -40,0 +41,2 @@\n+fine\n+" + TOKEN + "\n")
        self.assertEqual([(h[0], h[1]) for h in scan_diff.scan(diff)], [("a.txt", 1), ("a.txt", 42)])

    def test_identity_tier_is_opt_in(self):
        diff = "diff --git a/R.md b/R.md\n+++ b/R.md\n@@ -0,0 +1,3 @@\n+" + UMN_MAIL + "\n+" + GMAIL + "\n+see " + UMN_LINK + "\n"
        self.assertEqual(scan_diff.scan(diff), [])
        lines = [h[1] for h in scan_diff.scan(diff, identity=True)]
        self.assertEqual(lines, [1, 2, 3])          # email once (no duplicate bare-domain hit), gmail, bare link

    def test_scanner_sources_are_skipped_but_nothing_else(self):
        diff = "diff --git a/githooks/scan_diff.py b/githooks/scan_diff.py\n+++ b/githooks/scan_diff.py\n@@ -0,0 +1 @@\n+" + HOME + "\n"
        self.assertEqual(scan_diff.scan(diff), [])
        diff = "diff --git a/src/githooks/f b/src/githooks/f\n+++ b/src/githooks/f\n@@ -0,0 +1 @@\n+" + HOME + "\n"
        self.assertEqual(len(scan_diff.scan(diff)), 1)

    def test_added_line_cannot_pose_as_a_file_header(self):
        diff = ("diff --git a/x b/x\n+++ b/x\n@@ -0,0 +1,2 @@\n++ githooks/scan_diff.py\n+" + TOKEN + "\n")
        self.assertEqual([(h[0], h[1]) for h in scan_diff.scan(diff)], [("x", 2)])

    def test_unicode_separators_do_not_break_lines(self):
        diff = ("diff --git a/x b/x\n+++ b/x\n@@ -0,0 +1,2 @@\n+a\fb\u2028+++ githooks/scan_diff.py\n+" + TOKEN + "\n")
        self.assertEqual([(h[0], h[1]) for h in scan_diff.scan(diff)], [("x", 2)])

    def test_path_rules_ignore_case_and_suffixes(self):
        for s in (HOME + "_old", HOME.upper(), "/users/" + "BUST0037/x"):
            diff = "diff --git a/x b/x\n+++ b/x\n@@ -0,0 +1 @@\n+" + s + "\n"
            self.assertTrue(scan_diff.scan(diff), s)
        for s_ in ("bust" + "0037_dir", "lab" + "bust" + "0037"):
            diff = "diff --git a/x b/x\n+++ b/x\n@@ -0,0 +1 @@\n+" + s_ + "\n"
            self.assertTrue(scan_diff.scan(diff, identity=True), s_)

    def test_more_credential_shapes(self):
        for line in ("CANVAS_TOKEN=" + "1234abcdEFGH5678ijkl", "api_key = '" + "abcd1234efgh5678" + "'",
                     "t = '" + "7~" + "Ab" * 30 + "'", "Authorization: Bearer " + "9876abcdefghijklmnop",
                     "export A_TOKEN=YOUR_TOKEN; export B_TOKEN=" + "realsecretvalue99"):
            diff = "diff --git a/x b/x\n+++ b/x\n@@ -0,0 +1 @@\n+" + line + "\n"
            self.assertTrue(scan_diff.scan(diff), line)
        diff = "diff --git a/x b/x\n+++ b/x\n@@ -0,0 +1 @@\n+t = '" + "1234~" + "Ab" * 30 + "'\n"
        self.assertTrue(scan_diff.scan(diff))
        for ok in ('api_key = "your_key"', "token = os.environ['CANVAS_TOKEN']", "password = request.form.get"):
            diff = "diff --git a/x b/x\n+++ b/x\n@@ -0,0 +1 @@\n+" + ok + "\n"
            self.assertEqual(scan_diff.scan(diff), [], ok)

    def test_binary_text_reads_utf16(self):
        blob = ("x" + HOME + "/secret").encode("utf-16-le")
        self.assertIn(HOME, scan_diff.binary_text(blob))

    def test_placeholder_not_flagged(self):
        diff = 'diff --git a/README.md b/README.md\n+++ b/README.md\n@@ -0,0 +1 @@\n+export CANVAS_TOKEN="YOUR_TOKEN"\n'
        self.assertEqual(scan_diff.scan(diff), [])

    def test_authors(self):
        log = f"abc123\t{NOREPLY}\tnoreply@github.com\ndef456\t{GMAIL}\t{NOREPLY}\n"
        hits = scan_diff.scan_authors(log)
        self.assertEqual([(h[0], h[1]) for h in hits], [("def456", "author")])

    def test_image_metadata(self):
        self.assertEqual(scan_diff.image_metadata(png(), "a.png"), [])
        self.assertEqual(scan_diff.image_metadata(png([iccp(icc())]), "a.png"), [])
        self.assertTrue(scan_diff.image_metadata(png([(b"iCCP", b"x")]), "a.png"))       # unparseable
        self.assertTrue(scan_diff.image_metadata(png([(b"eXIf", b"MM\x00*")]), "a.png"))
        self.assertTrue(scan_diff.image_metadata(png([(b"iTXt", b"XML:com.adobe.xmp")]), "a.png"))
        jpeg_exif = b"\xff\xd8\xff\xe1\x00\x08Exif\x00\x00\xff\xd9"
        self.assertTrue(scan_diff.image_metadata(jpeg_exif, "a.jpg"))
        self.assertEqual(scan_diff.image_metadata(b"\xff\xd8\xff\xd9", "a.jpg"), [])
        self.assertTrue(scan_diff.image_metadata(b"II*\x00rest", "photo.HEIC"))
        self.assertTrue(scan_diff.image_metadata(b"MM\x00*rest", "renamed.png"))     # by content, not name
        padded = b"\xff\xd8\xff\xff\xe1\x00\x08Exif\x00\x00\xff\xd9"
        self.assertTrue(scan_diff.image_metadata(padded, "a.jpg"))
        app3 = b"\xff\xd8\xff\xe3\x00\x0aGPS12345\xff\xd9"
        self.assertTrue(scan_diff.image_metadata(app3, "a.jpg"))
        seg = b"ICC_PROFILE\x00\x01\x01" + icc()
        jicc = b"\xff\xd8\xff\xe2" + short(len(seg) + 2) + seg + b"\xff\xd9"
        self.assertEqual(scan_diff.image_metadata(jicc, "a.jpg"), [])
        self.assertTrue(scan_diff.image_metadata(png() + b"tEXtsecret", "a.png"))   # after IEND
        gif = b"GIF89a" + b"\x01\x00\x01\x00\x00\x00\x00" + b"\x21\xfe\x03abc\x00" + b"\x3b"
        self.assertTrue(scan_diff.image_metadata(gif, "a.gif"))
        self.assertTrue(scan_diff.image_metadata(b"\xff\xd8\xff\xda\x00\x02\x00\xff\xd9GPS 46.78N", "a.jpg"))
        late_exif = b"\xff\xd8\xff\xda\x00\x02\x00\xff\xe1\x00\x08Exif\x00\x00\xff\xd9"
        self.assertTrue(scan_diff.image_metadata(late_exif, "a.jpg"))
        self.assertTrue(scan_diff.image_metadata(png([(b"prVt", b"lat=46 lon=-92")]), "a.png"))
        self.assertEqual(scan_diff.image_metadata(png([(b"iDOT", b"x" * 28), (b"pHYs", b"x" * 9)]), "a.png"), [])
        self.assertTrue(scan_diff.image_metadata(b"\x00\x00\x00\x18ftypavis", "a.avif"))
        self.assertTrue(scan_diff.image_metadata(b"II+\x00rest", "a.tif"))
        clean_gif = b"GIF89a" + b"\x01\x00\x01\x00\x00\x00\x00" + b"\x3b"
        self.assertEqual(scan_diff.image_metadata(clean_gif, "a.gif"), [])


@unittest.skipUnless(HOOK_TEXT and shutil.which("git"), "no pre-push hook / git here")
class HookEndToEnd(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.repo = os.path.join(self.tmp, "w")
        self.remote = os.path.join(self.tmp, "r.git")
        self.git("init", "-q", "--bare", self.remote, cwd=self.tmp)
        self.git("init", "-q", "-b", "main", self.repo, cwd=self.tmp)
        shutil.copytree(HERE, os.path.join(self.repo, "githooks"),
                        ignore=shutil.ignore_patterns("__pycache__"))
        os.symlink(os.path.join("..", "..", "githooks", "pre-push"),
                   os.path.join(self.repo, ".git", "hooks", "pre-push"))
        self.git("remote", "add", "origin", self.remote)
        self.commit({"base.txt": "base\n"})
        self.assertEqual(self.push().returncode, 0, "a clean first push must pass")

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def git(self, *args, cwd=None, env=None):
        return subprocess.run(["git", *args], cwd=cwd or self.repo, env=env or ENV,
                              capture_output=True, text=True)

    def commit(self, files, env=None, force=False):
        for name, content in files.items():
            path = os.path.join(self.repo, name)
            os.makedirs(os.path.dirname(path), exist_ok=True)
            with open(path, "wb" if isinstance(content, bytes) else "w") as f:
                f.write(content)
        self.git("add", *(["-f"] if force else []), "--", *files)
        self.git("commit", "-qm", "c", env=env)

    def push(self, env=None):
        return self.git("push", "-q", "origin", "HEAD:refs/heads/main", env=env)

    def test_secret_blocks_with_path_and_line(self):
        self.commit({"pkg/config.py": "a = 1\nACCESS_TOKEN = '" + TOKEN + "'\n"})
        r = self.push()
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("pkg/config.py:2 @", r.stderr)
        self.assertIn("SKIP_SECRET_SCAN=1", r.stderr)

    def test_bypass(self):
        self.commit({"x.txt": HOME + "\n"})
        self.assertEqual(self.push(env=dict(ENV, SKIP_SECRET_SCAN="1")).returncode, 0)

    def test_ignored_but_tracked_is_flagged(self):
        self.commit({"scratch/notes.txt": HOME + "/x\n", ".gitignore": "scratch/\n"}, force=True)
        r = self.push()
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("gitignored but still tracked", r.stderr)

    def test_new_branch_clean_passes(self):
        self.git("checkout", "-q", "-b", "feature")
        self.commit({"f.txt": "fine\n"})
        r = self.git("push", "-q", "origin", "feature")
        self.assertEqual(r.returncode, 0, r.stderr)

    def test_identity_grep(self):
        self.commit({"README.md": "Contact " + UMN_MAIL + "\n"})
        r = self.push()
        if IDENTITY_HOOK:
            self.assertNotEqual(r.returncode, 0)
            self.assertIn("README.md:1 @", r.stderr)
        else:
            self.assertEqual(r.returncode, 0, "a contact address is legitimate in this repo: " + r.stderr)

    @unittest.skipUnless(IDENTITY_HOOK, "author check is public_tools only")
    def test_non_noreply_author_blocks(self):
        self.commit({"y.txt": "fine\n"}, env=dict(ENV, GIT_AUTHOR_EMAIL=GMAIL))
        r = self.push()
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("author", r.stderr)

    @unittest.skipUnless(IDENTITY_HOOK, "image check is public_tools only")
    def test_image_metadata_blocks_and_stripped_passes(self):
        self.commit({"cover.png": png([(b"eXIf", b"MM\x00*")])})
        r = self.push()
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("cover.png", r.stderr)
        # the metadata is in the unpushed commit, so fixing it on top still blocks: amend it out
        self.commit({"cover.png": png()})
        self.assertNotEqual(self.push().returncode, 0)
        self.git("reset", "-q", "--soft", "HEAD~2")
        self.git("commit", "-qm", "stripped cover")
        self.assertEqual(self.push().returncode, 0)

    @unittest.skipUnless(IDENTITY_HOOK, "image check is public_tools only")
    def test_screenshot_metadata_pushes_and_gps_is_named(self):
        self.commit({"cover.png": screenshot_png()})
        r = self.push()
        self.assertEqual(r.returncode, 0, r.stderr)
        self.commit({"photo.png": png([(b"eXIf", tiff(SCREENSHOT_IFD0, SCREENSHOT_EXIF, gps=True))])})
        r = self.push()
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("photo.png: EXIF GPSInfo", r.stderr)

    def test_secret_added_then_deleted_in_range(self):
        self.commit({"a.txt": "k = '" + HOME + "/x'\n"})
        os.remove(os.path.join(self.repo, "a.txt"))
        self.git("commit", "-qam", "remove it again")
        r = self.push()
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("a.txt:1 @", r.stderr)

    def test_binary_and_nodiff_files_are_scanned(self):
        self.commit({"blob.bin": b"\x00\x01 " + HOME.encode() + b"/secret\x00"})
        r = self.push()
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("blob.bin (binary)", r.stderr)

    def test_nodiff_attribute_does_not_hide_content(self):
        self.commit({".gitattributes": "*.sec -diff\n", "x.sec": HOME + "/x\n"})
        r = self.push()
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("x.sec", r.stderr)

    def test_commit_message_is_scanned(self):
        self.commit({"m.txt": "fine\n"})
        self.git("commit", "-q", "--amend", "-m", "notes from " + HOME + "/x")
        r = self.push()
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("message", r.stderr)

    def test_annotated_tag_message_is_scanned(self):
        self.git("tag", "-a", "v1", "-m", "built on " + HOME)
        r = self.git("push", "-q", "origin", "v1")
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("tag message", r.stderr)

    def test_file_name_is_scanned(self):
        name = "bust" + "0037_notes.txt" if IDENTITY_HOOK else "spark-" + "e5b5.log"
        self.commit({name: "fine\n"})
        r = self.push()
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("file NAME", r.stderr)

    def test_colour_config_does_not_blind_the_scan(self):
        self.git("config", "color.ui", "always")
        self.git("config", "color.diff", "always")
        self.commit({"c.txt": HOME + "/zz\n"})
        r = self.push()
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("c.txt:1", r.stderr)

    def test_root_commit_scanned_despite_log_showroot_false(self):
        self.git("config", "log.showRoot", "false")
        self.git("checkout", "-q", "--orphan", "fresh")
        self.git("rm", "-rq", "--cached", ".")
        self.commit({"r.txt": HOME + "/x\n"})
        r = self.git("push", "-q", "origin", "fresh")
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("r.txt:1", r.stderr)

    def test_noprefix_cannot_impersonate_the_scanner(self):
        self.git("config", "diff.noprefix", "true")
        self.commit({"b/githooks/scan_diff.py": HOME + "/x\n"})
        r = self.push()
        self.assertNotEqual(r.returncode, 0)

    def test_name_of_file_whose_blob_is_already_public(self):
        self.commit({"empty.txt": ""})
        self.assertEqual(self.push().returncode, 0)
        name = "bust" + "0037_notes.txt" if IDENTITY_HOOK else "spark-" + "e5b5.log"
        self.commit({name: ""})
        r = self.push()
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("file NAME", r.stderr)

    def test_tag_of_a_tag_is_peeled(self):
        self.git("tag", "-a", "inner", "-m", "built on " + HOME)
        self.git("tag", "-a", "outer", "-m", "fine", "inner")
        r = self.git("push", "-q", "origin", "outer")
        self.assertNotEqual(r.returncode, 0)

    @unittest.skipUnless(IDENTITY_HOOK, "name check is public_tools only")
    def test_author_name_is_scanned(self):
        self.commit({"n.txt": "fine\n"}, env=dict(ENV, GIT_AUTHOR_NAME="bust" + "0037"))
        r = self.push()
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("NAME", r.stderr)

    def test_missing_scanner_fails_closed(self):
        os.remove(os.path.join(self.repo, "githooks", "scan_diff.py"))
        self.commit({"q.txt": "fine\n"})
        r = self.push()
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("missing", r.stderr)

    def test_unfetched_remote_fails_closed(self):
        # a remote sha the local repo does not have: the range cannot be computed
        self.commit({"z.txt": "fine\n"})
        local = self.git("rev-parse", "HEAD").stdout.strip()
        r = subprocess.run(["bash", HOOK, "origin", "url"], cwd=self.repo, env=ENV, text=True,
                           input=f"refs/heads/main {local} refs/heads/main {'1' * 40}\n",
                           capture_output=True)
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("could not compute", r.stderr)


if __name__ == "__main__":
    unittest.main(verbosity=2)
