using BioSequences
using BioAlignments
using StatsBase
using FASTX

function ReadIn(fs::String)::Vector{LongSequence{DNAAlphabet{4}}}
    seq_lst = LongSequence{DNAAlphabet{4}}[]
    for line in eachline(fs)
        push!(seq_lst, LongSequence{DNAAlphabet{4}}(line))
    end
    seq_lst
end

function guina(seqLst::Vector{LongSequence{DNAAlphabet{4}}})::Dict{LongSequence{DNAAlphabet{4}},Int64}
    countmap(seqLst)
end

function global_align(query::LongSequence{DNAAlphabet{4}}, ref::LongSequence{DNAAlphabet{4}};gap_open=-10,gap_extend=-1)::PairwiseAlignment
    scoremodel = AffineGapScoreModel(EDNAFULL, gap_open=gap_open, gap_extend=gap_extend)
    res = pairalign(GlobalAlignment(), query, ref, scoremodel)
    aln = alignment(res)
    aln
end

function find_Deletions(aln::PairwiseAlignment, ref::LongSequence{DNAAlphabet{4}})::Vector{String}
    dl_len = zero(Int64)
    deletions = String[]
    first_p = 0
    for i in 1:length(ref)
        (_, op) = ref2seq(aln, i)
        if isdeleteop(op)
            dl_len += 1
            first_p = i
        else
            if dl_len != 0
                push!(deletions, "$(first_p - dl_len + 1)D$(dl_len):$(ref[(first_p - dl_len + 1):(first_p)])")
            end
            dl_len = 0
        end
    end
    if dl_len != 0
        push!(deletions, "$(first_p - dl_len + 1)D$(dl_len):$(ref[(first_p - dl_len + 1):(first_p)])")
    end
    deletions
end

function find_Insertions(aln::PairwiseAlignment, query::LongSequence{DNAAlphabet{4}})::Vector{String}
    ins_l = 0
    first_p = 0
    insertions = String[]
    k = 0
    for i in 1:length(query)
        (site, op) = seq2ref(aln, i)
        if isinsertop(op)
            ins_l += 1
            first_p = site
        else
            if ins_l != 0
                push!(insertions, "$(first_p)I$(ins_l):$(query[(i - ins_l):(i - 1)])")
            end
            ins_l = 0
        end
        k = i
    end
    if ins_l != 0
        push!(insertions, "$(first_p)I$(ins_l):$(query[(k - ins_l):(k - 1)])")
    end
    insertions
end

function find_Mismatch(aln::PairwiseAlignment, ref::LongSequence{DNAAlphabet{4}}, query::LongSequence{DNAAlphabet{4}})::Vector{String}
    mismatchs = []
    for i in 1:length(ref)
        (p, op) = ref2seq(aln, i)
        if op == Operation('X')
            push!(mismatchs, "$i:$(ref[i])2$(query[p])")
        end
    end
    mismatchs
end

function get_aligned_ref(aln::PairwiseAlignment)::LongSequence{DNAAlphabet{4}}
    LongSequence{DNAAlphabet{4}}([y for (_, y) in aln])
end

function get_aligned_query(aln::PairwiseAlignment)::LongSequence{DNAAlphabet{4}}
    LongSequence{DNAAlphabet{4}}([x for (x, _) in aln])
end

function calling(querys_mp::Dict, ref::LongSequence{DNAAlphabet{4}}, io::IO)
    println(io, "reference\tquery\tcount\tinsertions\tdeletions\tmismatches")
    for (query, v) in querys_mp
        aln = global_align(query, ref)
        ins_lst = join(find_Insertions(aln, query), "|")
        del_lst = join(find_Deletions(aln, ref), "|")
        mismatchs = join(find_Mismatch(aln, ref, query), "|")
        rf = get_aligned_ref(aln)
        qy = get_aligned_query(aln)
        println(io, "$(rf)\t$(qy)\t$(v)\t$(ins_lst)\t$(del_lst)\t$(mismatchs)")
    end
end

function main(;reference::String="none", inPath::String="none",outPath::String="none",start::Int64=118, fin::Int64=155)
    refs = open(FASTA.Reader, reference, index = "$reference.fai")
    all_lib = collect(keys(refs.index.names))
    if isdir(outPath) == false
        mkpath(outPath)
    end
    for i in all_lib
        if isfile("$inPath/$i.filter.out") && stat("$inPath/$i.filter.out").size > 0
            open("$outPath/$i.indelProfile.tsv", "w+") do IO
                calling(guina(ReadIn("$inPath/$i.filter.out")), LongSequence{DNAAlphabet{4}}(FASTX.sequence(refs[i])[start:fin]), IO)
            end
        else
            println("$inPath/$i.filter.out is empty!")
        end
    end
end

main(reference=snakemake.config["ref"],inPath=snakemake.input[1],outPath=snakemake.output[1],start=snakemake.config["RealignToTrap"]["start"],fin=snakemake.config["RealignToTrap"]["end"])