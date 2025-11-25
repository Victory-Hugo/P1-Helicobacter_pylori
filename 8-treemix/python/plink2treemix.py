#!/usr/bin/python3
import sys, gzip

if len(sys.argv) < 3:
    print("plink2treemix.py [gzipped input file] [gzipped output file]")
    sys.exit(1)

infile  = gzip.open(sys.argv[1], "rt", encoding="utf-8")
outfile = gzip.open(sys.argv[2], "wt", encoding="utf-8")

pop2rs = {}
rss_order = []            # 保持 SNP 输出顺序
seen_rs  = set()

# 跳过前两行
infile.readline()
infile.readline()

for line in infile:
    cols = line.strip().split()
    if len(cols) < 8:
        continue
    chr_id, rs, pop = cols[0], cols[1], cols[2]
    mc, total = map(int, cols[6:8])

    snp_id = f"{chr_id}_{rs}"
    if snp_id not in seen_rs:
        rss_order.append(snp_id)
        seen_rs.add(snp_id)

    pop2rs.setdefault(pop, {})
    if snp_id in pop2rs[pop]:
        # 如遇重复，同群体同 SNP 累加
        old_mc, old_total = map(int, pop2rs[pop][snp_id].split())
        mc    += old_mc
        total += old_total
    pop2rs[pop][snp_id] = f"{mc} {total}"

# 写文件
pops = list(pop2rs.keys())
print(*pops, file=outfile)

for snp_id in rss_order:
    row = []
    for pop in pops:
        if snp_id in pop2rs[pop]:
            mc, total = map(int, pop2rs[pop][snp_id].split())
            row.append(f"{mc},{total-mc}")
        else:
            row.append("0,0")
    print(" ".join(row), file=outfile)
