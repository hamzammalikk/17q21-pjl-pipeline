import json, subprocess, urllib.request
out = subprocess.run(["bcftools", "query", "-f", "%CHROM\t%POS\t%REF\t%ALT\n", "results/cand.sites.vcf.gz"], capture_output=True, text=True, check=True).stdout.split("\n")
variants = []
for line in out:
    if not line.strip():
        continue
    c, p, ref, alt = line.split("\t")
    variants.append(f"{c.replace('chr', '')} {p} . {ref} {alt} . . .")
print("sending", len(variants), "variants")
req = urllib.request.Request(
    "https://rest.ensembl.org/vep/human/region?canonical=1&regulatory=1",
    data=json.dumps({"variants": variants}).encode(),
    headers={"Content-Type": "application/json", "Accept": "application/json"})
res = json.load(urllib.request.urlopen(req, timeout=600))
with open("results/vep.tsv", "w") as f:
    f.write("Location\tAlleles\tExisting_variation\tConsequence\tSYMBOL\tIMPACT\tREGULATORY\n")
    for r in res:
        allt = r.get("transcript_consequences", [])
        tc = [t for t in allt if t.get("canonical") == 1] or allt
        genes = ",".join(sorted({t["gene_symbol"] for t in tc if t.get("gene_symbol")})) or "-"
        ms = r["most_severe_consequence"]
        imp = next((t.get("impact") for t in allt if ms in t.get("consequence_terms", [])), "-")
        rs = ",".join(c["id"] for c in r.get("colocated_variants", []) if str(c.get("id", "")).startswith("rs")) or "-"
        reg = ",".join(sorted({t for x in r.get("regulatory_feature_consequences", []) for t in x.get("consequence_terms", [])})) or "-"
        f.write(f"{r['seq_region_name']}:{r['start']}\t{r['allele_string']}\t{rs}\t{ms}\t{genes}\t{imp}\t{reg}\n")
print("wrote results/vep.tsv with", len(res), "rows")
