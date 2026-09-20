import json, subprocess, urllib.request, collections
q = subprocess.run(["bcftools", "query", "-f", "%CHROM\t%POS\t%REF\t%ALT\n", "results/cand.sites.vcf.gz"], capture_output=True, text=True, check=True).stdout
keys = set()
for line in q.strip().split("\n"):
    c, p, ref, alt = line.split("\t")
    keys.add(f"{c}_{p}_{ref}_{alt}_b38")
print("candidates:", len(keys))
hits = []
for tis in ["Lung", "Whole_Blood"]:
    url = "https://gtexportal.org/api/v2/association/singleTissueEqtlByLocation?tissueSiteDetailId=" + tis + "&start=39700000&end=40000000&chromosome=chr17&datasetId=gtex_v8"
    d = json.load(urllib.request.urlopen(url, timeout=300))["singleTissueEqtl"]
    print(tis, "significant eQTL rows in the whole window:", len(d))
    hits += [r for r in d if r["variantId"] in keys]
hits.sort(key=lambda r: (r["pos"], r["tissueSiteDetailId"], r["geneSymbol"]))
with open("results/eqtl_hits.tsv", "w") as f:
    f.write("pos\tsnpId\tvariantId\ttissue\tgene\tnes\tpValue\n")
    for r in hits:
        f.write(f"{r['pos']}\t{r['snpId']}\t{r['variantId']}\t{r['tissueSiteDetailId']}\t{r['geneSymbol']}\t{r['nes']}\t{r['pValue']}\n")
print("candidates that are eQTLs:", len({r["variantId"] for r in hits}), "of", len(keys))
for k, v in sorted(collections.Counter((r["tissueSiteDetailId"], r["geneSymbol"]) for r in hits).items()):
    print(k, v)
print("--- the 4 protein-altering variants and the lead SNP ---")
focus = {39872381, 39905943, 39905964, 39908216, 39913696}
for r in hits:
    if r["pos"] in focus:
        print(r["pos"], r["snpId"], r["tissueSiteDetailId"], r["geneSymbol"], r["nes"], r["pValue"])
