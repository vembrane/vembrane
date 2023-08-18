import pysam

with pysam.VariantFile("test.bcf", "r") as f:
    existing = set(f.header.info.keys())

    for variant in f:
        existing -= set(variant.info.keys())

    for e in existing:
        print(f"##INFO=<ID={e}")

# with pysam.VariantFile("test.bcf", "r") as f:

# for e in header.info:
#     header.info.clear_header()

# print([x for x in header.info])

# with pysam.VariantFile("test.clean.bcf", "o", header=f.header) as o:
# existing = set(f.header.info.keys())

# for variant in f:
#     existing -= set(variant.info.keys())
