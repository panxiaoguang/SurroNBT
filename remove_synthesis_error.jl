using Fire

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

function removeBAC(ref::String, ctrl::Vector, case::Vector, start::Int64)
    filter(x -> (!((x in ctrl) && x != ref[start:end])), case)
end

@main function fin(;reference::String="none", outPath::String="none", start::Int64=114)
    refs = readRef(reference)
    all_lib = collect(keys(refs))
    for i in all_lib
        if isfile("$(spCas9)/$i.out") && isfile("$(WT)/$i.out")
            case = ReadIn("$(spCas9)/$i.out")
            ctrl = ReadIn("$(WT)/$i.out")
            open("$outPath/$i.filter.out", "w+") do IO
                for line in removeBAC(refs[i], ctrl, case, start)
                    write(IO, line, "\n")
                end
            end
        else
            println("$i have been removed!")
        end
    end
end

