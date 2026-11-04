import mudata as md

path = "/data/20_data_ayat/scenicplus_workflow/scenicplus_out/WT/ACC_GEX.h5mu"

mdata = md.read_h5mu(path)

rna = mdata.mod["scRNA"]

print("before unique:", rna.var_names.is_unique)

#Forcer d'enlever les doublons, sinon SCENIC+ bloque
rna = rna[:, ~rna.var_names.duplicated()].copy()
rna.var_names_make_unique()

print("after unique:", rna.var_names.is_unique)

mdata.mod["scRNA"] = rna

mdata.write(path)
import mudata as md

m = md.read_h5mu("/data/20_data_ayat/scenicplus_workflow/scenicplus_out/WT/ACC_GEX.h5mu")

rna = m.mod["scRNA"]

print("unique:", rna.var_names.is_unique)
print("duplicates:", sum(rna.var_names.duplicated()))
print("example dup:", rna.var_names[rna.var_names.duplicated()][:10])
