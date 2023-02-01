using BioSequences
using BioAlignments
using StatsBase
using Fire

function ReadIn(fs::String)::Vector{LongDNASeq}
    seq_lst = LongDNASeq[]
    for line in eachline(fs)
        push!(seq_lst, LongDNASeq(line))
    end
    seq_lst
end

function readRef(fs::String)::Dict{String,LongDNASeq}
    seq_name = ""
    fasta = Dict{String,LongDNASeq}()
    for line in eachline(fs)
        if startswith(line, ">")
            seq_name = replace(line, ">" => "")
        else
            fasta[seq_name] = LongDNASeq(line)
        end
    end
    fasta
end

function guina(seqLst::Vector{LongDNASeq})::Dict{LongSequence{DNAAlphabet{4}},Int64}
    countmap(seqLst)
end

function global_align(query::LongDNASeq, ref::LongDNASeq;gap_open=-5,gap_extend=-1)::PairwiseAlignment
    scoremodel = AffineGapScoreModel(EDNAFULL, gap_open=gap_open, gap_extend=gap_extend)
    res = pairalign(GlobalAlignment(), query, ref, scoremodel)
    aln = alignment(res)
    aln
end

function find_Deletions(aln::PairwiseAlignment, ref::LongDNASeq)::Vector{String}
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
                push!(deletions, "$(first_p - dl_len)D$(dl_len)")
            end
            dl_len = 0
        end
    end
    if dl_len != 0
        push!(deletions, "$(first_p - dl_len)D$(dl_len)")
    end
    deletions
end

function find_Insertions(aln::PairwiseAlignment, query::LongDNASeq)::Vector{String}
    ins_l = 0
    first_p = 0
    insertions = String[]
    for i in 1:length(query)
        (site, op) = seq2ref(aln, i)
        if isinsertop(op)
            ins_l += 1
            first_p = site
        else
            if ins_l != 0
                push!(insertions, "$(first_p)I$(ins_l)")
            end
            ins_l = 0
        end
    end
    if ins_l != 0
        push!(insertions, "$(first_p)I$(ins_l)")
    end
    insertions
end

function find_Mismatch(aln::PairwiseAlignment, ref::LongDNASeq)::Vector{String}
    mismatchs = []
    for i in 1:length(ref)
        (_, op) = ref2seq(aln, i)
        if op == Operation('X')
            push!(mismatchs, "$(i)M1")
        end
    end
    mismatchs
end

function get_aligned_ref(aln::PairwiseAlignment)::LongDNASeq
    LongDNASeq([y for (_, y) in aln])
end

function get_aligned_query(aln::PairwiseAlignment)::LongDNASeq
    LongDNASeq([x for (x, _) in aln])
end

function calling(querys_mp::Dict, ref::LongDNASeq, io::IO)
    println(io, "reference\tquery\tcount\tinsertions\tdeletions\tmismatches")
    for (query, v) in querys_mp
        aln = global_align(query, ref)
        ins_lst = join(find_Insertions(aln, query), "|")
        del_lst = join(find_Deletions(aln, ref), "|")
        mismatchs = join(find_Mismatch(aln, ref), "|")
        rf = get_aligned_ref(aln)
        qy = get_aligned_query(aln)
        println(io, "$(rf)\t$(qy)\t$(v)\t$(ins_lst)\t$(del_lst)\t$(mismatchs)")
    end
end

@main function main(;reference::String="none", inPath::String="none",outPath::String="none",start::Int64=114)
    refs = readRef(reference)
    all_lib = collect(keys(refs))
    for i in all_lib
        if isfile("$inPath/$i.filter.out")
            open("$outPath/$i.indelProfile.tsv", "w+") do IO
	            println("call indel for ", i)
                calling(guina(ReadIn("$inPath/$i.filter.out")), refs[i][start:end], IO)
            end
        end
    end
end


