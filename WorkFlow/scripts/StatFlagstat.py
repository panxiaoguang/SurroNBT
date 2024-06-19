def mapping_stat(flagstat):
    ## parse samtools flagstat and got mapping ratio
    f = open(flagstat, "r")
    data = f.readlines()
    total = int(data[0].split(" ")[0])
    mapped = int(data[6].split(" ")[0])
    return round(mapped / total * 100, 2)


if __name__ == "__main__":
    map_case = mapping_stat(snakemake.input["Case"])
    map_control = mapping_stat(snakemake.input["Control"])
    with open(snakemake.output[0], "w") as f:
        f.write("map_case\tmap_control\n")
        f.write(f"{map_case}\t{map_control}\n")
