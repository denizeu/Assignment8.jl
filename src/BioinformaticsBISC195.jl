module BioinformaticsBISC195

using Base: stop_reading
export normalizeDNA
export composition
export gc_content
export complement
export reverse_complement
export parse_fasta
# # uncomment the following line if you intend to use BioSequences types
# using BioSequences

"""
    normalizeDNA(::AbstractString)

Ensures that a sequence only contains valid bases
(or `'N'` for unknown bases).
Returns a String.
"""
function normalizeDNA(sequence)
    sequence= uppercase(sequence)
    rep = Dict('R' => 'N', 'Y' => 'N', 'S' => 'N', 'K' => 'N', 'M' => 'N', 'B' => 'N', 'D' => 'N', 'H' => 'N', 'V' => 'N', '.' => 'N', '-' => 'N', 'A' => 'A', 'G' => 'G', 'C' => 'C', 'T' => 'T', 'N' => 'N')
    return join([rep[c] for c in sequence])
end
    #=sequence = uppercase(string(sequence))
    for base in sequence # note: `N` indicates an unknown base
        occursin(base, "AGCTN") || error("invalid base $base")
    end
    return sequence # change to `return LongDNASeq(seq)` if you want to try to use BioSequences types
end=#


# Your code here.
#=function composition(sequence)
    sequence = normalizeDNA(sequence) #make uppercase string, check invalid bases
    a = c = g = t = 0 #sets all 4 variables to `0`
    comps = Dict() #initialize empty Dictionary for bases
    stopindex= length(sequence) #stops at length of sequence

    for i in 1:stopindex #for loop: traverses each sequence?
        for base in sequence  #tells what to look at: but does it recognize the bases?
        comp = sequence[i:stopindex] #indexing sequence?
        comp = normalizeDNA(comp)
        println(comp)
        if base == 'A'
            comps['A'] = a = a + 1
            println(comps['A'])
        end
        if base == 'G'
            g = g + 1
           end
        if base == 'T'
            t = t + 1 
        end
        if base == 'C'
            c = c + 1 
        end
    end
end=#
function composition(sequence)
    sequence = normalizeDNA(sequence) #make uppercase string, check invalid bases
    #a = c = g = t = n= 0 #sets all 4 variables to `0`
    comp = Dict()
    for i in 1:length(sequence)
        pos = sequence[i]
        if haskey(comp, pos) == true
            comp[pos] = comp[pos] + 1
        else
            comp[pos] = 1
        end
    end
    return comp
end
# Don't forget to export your functions!

function gc_content(sequence)
    ng = 0
    nc = 0
    sequence = normalizeDNA(sequence)
    for i in 1:length(sequence)
        if sequence[i] == 'G'
            ng = ng + 1
        elseif sequence[i] == 'C'
            nc = nc + 1
        end
    end
    return (ng + nc) / length(sequence)
end

function complement(sequence)
    sequence = normalizeDNA(sequence)
    ret = ""
    for i in 1:length(sequence)
        if sequence[i] == 'A'
            ret = ret * 'T'
        elseif sequence[i] == 'T'
            ret = ret * 'A'
        elseif sequence[i] == 'G'
            ret = ret * 'C'
        elseif sequence[i] == 'C'
            ret = ret * 'G'
        elseif sequence[i] == 'N'
            ret = ret * 'N'
        end
    end
    return ret 
end

function reverse_complement(sequence)
    st= ""
    sequence =reverse(sequence)
    for i in 1:length(sequence)
        st = st*string(complement(sequence[i]))
    end
    return st
end

function parse_fasta(path)
    headers= String[]
    sequences= []
    tmp= ""
    for line in eachline(path)
        if startswith(line, '>')
            if !isempty(tmp)
                push!(sequences, tmp)
            end
            tmp= ""
            push!(headers, SubString(line, 2))
        else
            line = normalizeDNA(line)
            tmp= tmp*line
        end
    end
        push!(sequences, tmp)
    return (headers, sequences)
end

