#!/usr/bin/env python3
"""Extract a named character vector from an R .rds file without R.

Why this exists: the authors deposit their hand curation of Lactobacillus V4
sequences as LactoSpecies_ByHand.rds -- a named character vector whose names
are species labels and whose values are 235 bp sequences. We need it as a TSV
before an R environment exists on the cluster, so this parses R's XDR
serialization directly. Only the SEXP types that occur in that file are
implemented; anything else raises rather than guessing.

Usage:
    decode_rds_strings.py <in.rds> <out.tsv>
    decode_rds_strings.py --selftest
"""
import gzip
import struct
import sys

NILVALUE_SXP = 254
REFSXP       = 255
SYMSXP       = 1
LISTSXP      = 2
CHARSXP      = 9
INTSXP       = 13
REALSXP      = 14
STRSXP       = 16
VECSXP       = 19


class RObj:
    """A parsed SEXP: its Python value plus any attributes (names, class...)."""
    def __init__(self, value, attr=None):
        self.value = value
        self.attr = attr or {}

    def __repr__(self):
        return "RObj(%r, attr=%r)" % (self.value, list(self.attr))


class Reader:
    def __init__(self, buf):
        self.b = buf
        self.i = 0
        self.refs = []          # symbol reference table, 1-indexed in R

    def int(self):
        v = struct.unpack(">i", self.b[self.i:self.i + 4])[0]
        self.i += 4
        return v

    def dbl(self):
        v = struct.unpack(">d", self.b[self.i:self.i + 8])[0]
        self.i += 8
        return v

    def raw(self, n):
        v = self.b[self.i:self.i + n]
        self.i += n
        return v

    def item(self):
        flags = self.int()
        typ = flags & 0xFF
        has_attr = bool(flags & (1 << 9))
        has_tag = bool(flags & (1 << 10))

        if typ == NILVALUE_SXP:
            return None

        if typ == REFSXP:
            idx = flags >> 8
            if idx == 0:
                idx = self.int()
            return self.refs[idx - 1]

        if typ == SYMSXP:
            name = self.item()
            sym = name.value if isinstance(name, RObj) else name
            self.refs.append(sym)
            return sym

        if typ == CHARSXP:
            n = self.int()
            if n == -1:
                return RObj(None)
            return RObj(self.raw(n).decode("utf-8", "replace"))

        if typ == LISTSXP:
            cell_attr = self.item() if has_attr else None
            tag = self.item() if has_tag else None
            car = self.item()
            cdr = self.item()
            return ("pairlist", tag, car, cdr, cell_attr)

        if typ == STRSXP:
            n = self.int()
            vals = []
            for _ in range(n):
                el = self.item()
                vals.append(el.value if isinstance(el, RObj) else el)
            obj = RObj(vals)
        elif typ == VECSXP:
            n = self.int()
            obj = RObj([self.item() for _ in range(n)])
        elif typ == INTSXP:
            n = self.int()
            obj = RObj([self.int() for _ in range(n)])
        elif typ == REALSXP:
            n = self.int()
            obj = RObj([self.dbl() for _ in range(n)])
        else:
            raise NotImplementedError(
                "SEXP type %d not implemented (byte offset %d)" % (typ, self.i))

        if has_attr:
            obj.attr = self._pairlist_to_dict(self.item())
        return obj

    def _pairlist_to_dict(self, node):
        out = {}
        while isinstance(node, tuple) and node[0] == "pairlist":
            _, tag, car, cdr, _ = node
            if tag is not None:
                out[tag] = car
            node = cdr
        return out


def read_rds(path):
    with open(path, "rb") as f:
        head = f.read(2)
    opener = gzip.open if head == b"\x1f\x8b" else open
    with opener(path, "rb") as f:
        buf = f.read()

    if buf[:2] != b"X\n":
        raise ValueError("not XDR-serialized R data (magic %r); "
                         "was the file saved with ascii=TRUE?" % buf[:2])
    r = Reader(buf)
    r.i = 2
    version = r.int()
    r.int()
    r.int()
    if version >= 3:
        enc_len = r.int()
        r.raw(enc_len)
    return r.item()


def named_character(obj):
    """Return [(name, value)] for a named character vector."""
    if not isinstance(obj, RObj) or not isinstance(obj.value, list):
        raise ValueError("top-level object is not a vector: %r" % (obj,))
    names_obj = obj.attr.get("names")
    if names_obj is None:
        raise ValueError("vector has no names attribute; attributes present: %s"
                         % sorted(obj.attr))
    names = names_obj.value
    if len(names) != len(obj.value):
        raise ValueError("names (%d) and values (%d) differ in length"
                         % (len(names), len(obj.value)))
    return list(zip(names, obj.value))


def _selftest():
    def i32(v):
        return struct.pack(">i", v)

    def charsxp(s):
        bs = s.encode()
        return i32(0x00040009) + i32(len(bs)) + bs

    def strsxp(items, has_attr=False):
        flags = 16 | ((1 << 9) if has_attr else 0)
        return i32(flags) + i32(len(items)) + b"".join(charsxp(x) for x in items)

    def symsxp(name):
        return i32(1) + charsxp(name)

    values = ["ACGT" * 3, "TTTT", "GGGGA"]
    names = ["Lactobacillus crispatus", "Lactobacillus iners", "NA_like"]

    payload = (
        b"X\n" + i32(2) + i32(0x040401) + i32(0x020300)
        + strsxp(values, has_attr=True)
        + i32(LISTSXP | (1 << 10)) + symsxp("names") + strsxp(names)
        + i32(NILVALUE_SXP)
    )

    import tempfile, os
    fd, path = tempfile.mkstemp(suffix=".rds")
    os.close(fd)
    with open(path, "wb") as f:
        f.write(payload)
    try:
        pairs = named_character(read_rds(path))
    finally:
        os.unlink(path)

    assert pairs == list(zip(names, values)), pairs
    print("selftest OK: parsed %d named entries" % len(pairs))


def main():
    if len(sys.argv) == 2 and sys.argv[1] == "--selftest":
        _selftest()
        return
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    src, dst = sys.argv[1], sys.argv[2]
    pairs = named_character(read_rds(src))

    # Which side holds the sequence is not knowable a priori: the authors key
    # this vector BY sequence (names) because the rownames of their taxonomy
    # table are sequences -- see tax[names(spc.lac.byhand), "Species"] in
    # RepRefine_PNAS.Rmd. Detect rather than assume, and refuse if ambiguous.
    def is_dna(s):
        return isinstance(s, str) and len(s) > 50 and not (set(s) - set("ACGTN"))

    n_name_dna = sum(1 for n, _ in pairs if is_dna(n))
    n_val_dna = sum(1 for _, v in pairs if is_dna(v))

    if n_name_dna == len(pairs) and n_val_dna == 0:
        rows = [(v, n) for n, v in pairs]       # names = sequence, values = species
        orientation = "names are sequences, values are species labels"
    elif n_val_dna == len(pairs) and n_name_dna == 0:
        rows = list(pairs)                      # names = species, values = sequence
        orientation = "names are species labels, values are sequences"
    else:
        sys.exit("cannot tell which side is the sequence: %d/%d names look like "
                 "DNA, %d/%d values do" % (n_name_dna, len(pairs),
                                           n_val_dna, len(pairs)))

    lengths = sorted({len(seq) for _, seq in rows})
    with open(dst, "w") as out:
        out.write("species\tlength\tsequence\n")
        for species, seq in rows:
            out.write("%s\t%d\t%s\n" % (species, len(seq), seq))

    print("%s -> %s" % (src, dst))
    print("  orientation    : %s" % orientation)
    print("  entries        : %d" % len(rows))
    print("  distinct labels: %d" % len({s for s, _ in rows}))
    print("  distinct seqs  : %d" % len({q for _, q in rows}))
    print("  sequence length: %s" % ", ".join(str(x) for x in lengths))


if __name__ == "__main__":
    main()
