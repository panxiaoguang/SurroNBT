function ReadIn(fs::String)
    seq_lst = String[]
    for line in eachline(fs)
        push!(seq_lst, line)
    end
    seq_lst
end

function readRef(fs::String)
    seq_name = ""
    fasta = Dict{String,String}()
    for line in eachline(fs)
        if startswith(line, ">")
            seq_name = replace(line, ">" => "")
        else
            fasta[seq_name] = line
        end
    end
    fasta
end

function removeBAC(ref::String, ctrl::Vector, case::Vector, start::Int64, fin::Int64)
    filter(x -> (!((x in ctrl) && x != ref[start:fin])), case)
end

function fin(;spCas9::String="none",WT::String="none",reference::String="none", outPath::String="none", start::Int64=118, fin::Int64=155)
    refs = readRef(reference)
    all_lib = collect(keys(refs))
    if isdir(outPath) == false
        mkdir(outPath)
    end
    for i in all_lib
        if isfile("$(spCas9)/$i.out") && isfile("$(WT)/$i.out")
            case = ReadIn("$(spCas9)/$i.out")
            ctrl = ReadIn("$(WT)/$i.out")
            open("$outPath/$i.filter.out", "w+") do IO
                for line in removeBAC(refs[i], ctrl, case, start, fin)
                    write(IO, line, "\n")
                end
            end
        else
            println("$i have been removed!")
        end
    end
end
fin(spCas9=snakemake.input["Case"], WT=snakemake.input["Control"], reference=snakemake.config["ref"], outPath=snakemake.output[1], start=snakemake.config["RemoveSyntheticErrors"]["start"],fin=snakemake.config["RemoveSyntheticErrors"]["end"])