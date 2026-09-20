import pandas as pd
ld = pd.read_csv("results/ld_lead.ld", sep=r"\s+")[["BP_B", "SNP_B", "R2"]]
ld.columns = ["pos", "plink_id", "r2_with_lead"]
ld = ld[ld.r2_with_lead >= 0.6]
vep = pd.read_csv("results/vep.tsv", sep="\t")
vep["plink_id"] = "17:" + vep["Location"].str.split(":").str[1] + ":" + vep["Alleles"].str.replace("/", ":")
vep = vep.drop_duplicates("plink_id")
af = pd.read_csv("results/pjl17.qc.afreq", sep="\t")[["ID", "ALT_FREQS"]]
af.columns = ["plink_id", "PJL_alt_freq"]
eq = pd.read_csv("results/eqtl_hits.tsv", sep="\t")
eq["plink_id"] = eq["variantId"].str.replace("chr", "").str.replace("_b38", "").str.replace("_", ":")
eq["tag"] = eq["tissue"].str.replace("Whole_Blood", "Blood") + ":" + eq["gene"] + eq["nes"].apply(lambda x: "+" if x > 0 else "-")
g = eq.groupby("plink_id")
eqs = pd.DataFrame({"eQTL": g["tag"].apply(lambda t: ";".join(sorted(t))), "n_eQTL": g.size(), "min_eQTL_p": g["pValue"].min()}).reset_index()
final = ld.merge(vep[["plink_id", "Existing_variation", "Consequence", "IMPACT", "SYMBOL", "REGULATORY"]], on="plink_id", how="left").merge(af, on="plink_id", how="left").merge(eqs, on="plink_id", how="left")
order = {"HIGH": 0, "MODERATE": 1, "LOW": 2, "MODIFIER": 3}
final["rank"] = final["IMPACT"].map(order)
final = final.sort_values(["rank", "r2_with_lead"], ascending=[True, False]).drop(columns="rank")
final["rsID"] = final["Existing_variation"].str.split(",").str[0]
final.to_csv("results/final_table.tsv", sep="\t", index=False)
print(len(final), "rows saved to results/final_table.tsv")
print(final[["pos", "rsID", "r2_with_lead", "PJL_alt_freq", "Consequence", "SYMBOL", "n_eQTL"]].head(12).to_string(index=False))
