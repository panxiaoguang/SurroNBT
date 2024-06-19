using CSV
using DataFrames

function panduan(x::Vector{T},buyao::Vector{T}) where T <: AbstractString
    !all(in.(x, Ref(buyao)))
end

function modify_df(fs::String,lftS::Int64,lftE::Int64,rgtS::Int64,rgtE::Int64)
    df = DataFrame(CSV.File(fs, delim="\t"))
    buyao = vcat(string.(lftS:lftE) .* "I1", string.(rgtS:rgtE) .* "I1", string.(lftS:lftE) .* "D1", string.(rgtS:rgtE) .* "D1")
    if (eltype(df.insertions) == Missing && eltype(df.deletions) == Missing) | (size(df)[1] == 0)
        return Missing
    elseif eltype(df.insertions) == Missing
        df.insertions = string.(df.insertions)
        df.insertions .= ""
    elseif eltype(df.deletions) == Missing
        df.deletions = string.(df.deletions)
        df.deletions .= ""
    end
    df[:,[:insertions,:deletions]] .= ifelse.(ismissing.(df[:,[:insertions,:deletions]]), "", df[:,[:insertions,:deletions]])
    df = select(df, :reference, :query, :count, AsTable(4:5) => ByRow(x -> replace.(split(strip(join([x.insertions,x.deletions], "|"), ['|']), "|"),r":[AGCTN]+" => "")) => :indels)
    sort!(df, order(:count, rev=true))
    df = df[panduan.(df.indels,Ref(buyao)),:]
    df
end

buxing(x) = !((length(x) == 1) && (x[1] == ""))

function calIndel(x)
    if x == Missing
        return 0, 0
    end
    cleanReadNumbers = sum(x.count)
    indelReadNumbers = sum(x[buxing.(x.indels),:].count)
    cleanReadNumbers, indelReadNumbers
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

function main(;reference::String="none",inPath::String="none",out::String="none",lftS::Int64=1,lftE::Int64=25,rgtS::Int64=32,rgtE::Int64=37)
    outpath = dirname(out)
    if !isdir(outpath)
        mkdir(outpath)
    end
    f = open(out, "w")
    refs = readRef(reference)
    all_lib = collect(keys(refs))
    println(f,"Label\tCleanReadNumbers\tIndelReadNumbers")
    for i in all_lib
        if isfile("$inPath/$i.indelProfile.tsv")
            clean, indels = calIndel(modify_df("$inPath/$i.indelProfile.tsv",lftS,lftE,rgtS,rgtE))
            println(f, i, "\t", clean, "\t", indels)
        end
    end
    close(f)
end

main(reference=snakemake.config["ref"],inPath=snakemake.input[1],out=snakemake.output[1],lftS=snakemake.config["CreateIndelCountFile"]["leftStart"],lftE=snakemake.config["CreateIndelCountFile"]["leftEnd"],rgtS=snakemake.config["CreateIndelCountFile"]["rightStart"],rgtE=snakemake.config["CreateIndelCountFile"]["rightEnd"])