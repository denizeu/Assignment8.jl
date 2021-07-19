module BioinformaticsBISC195

using Base: stop_reading
export normalizeDNA
export composition
export gc_content
export complement
export reverse_complement
export parse_fasta
export uniqueKmers
export kmerdist

# ### 1. NormalizeDNA Function
"""
    normalizeDNA(::AbstractString)

Ensures that a sequence only contains valid bases
(or `'N'` for unknown bases).
Returns a String.

Examples  
≡≡≡≡≡≡≡≡≡≡

    julia> normalizeDNA("acatg")
    "ACATG"

    julia> normalizeDNA("LAC")
    ERROR: key 'L' not found

    julia> normalizeDNA(7)
    ERROR: MethodError

    julia> normalizeDNA('G')
    "G"

    julia> s = normalizeDNA('C');

    julia> typeof(s)
    String
"""
function normalizeDNA(sequence)
    sequence= uppercase(sequence)
    rep = Dict('R' => 'N', 'Y' => 'N', 'S' => 'N', 'K' => 'N', 'M' => 'N', 'B' => 'N', 'D' => 'N', 'H' => 'N', 'W' => 'N','V' => 'N', '.' => 'N', '-' => 'N', 'A' => 'A', 'G' => 'G', 'C' => 'C', 'T' => 'T', 'N' => 'N')
    return join([rep[c] for c in sequence])
end

# ### 2. Composition Function
"""
    composition(sequence)

Counts the number of each type of base
in a DNA sequence and returns a Dictionary of how many of each exist
in the order A, C, G, T.
Converts any ambiguous bases not in AGCT to "N".

Examples  
≡≡≡≡≡≡≡≡≡≡

    julia> composition("ACGGG")
    'A' => 1
    'G' => 3
    'C' => 1

    julia> composition("ACC")
    'A' => 1
    'C' => 2

    julia> A,C,T = composition("acctt")
    'A' => 1
    'T' => 2
    'C' => 2

    julia> composition("BAGGGRRRR")
    'A' => 1
    'G' => 3
    'N' => 5   
"""
function composition(sequence)
    sequence = normalizeDNA(sequence) #make uppercase string, check invalid bases
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

# ### 3. GC_Content Function
"""
    gc_content(sequence)

Calculates the GC ratio of a DNA sequence.
The GC ratio is the total number of G and C bases divided by the total length of the sequence.

Examples  
≡≡≡≡≡≡≡≡≡≡

    julia> gc_content("CATA")
    0.25

    julia> gc_content("cccg") * 20
    20.0

    julia> gc_content("AaTt")
    0.0

    julia> gc_content("ATty")
    0.0
"""
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

# ###4. Complement Function
"""
    complement(base)

Get the DNA complement of the provided base:

    A <-> T
    G <-> C

Accepts uppercase or lowercase `String` or `Char`,
but always returns an uppercase `String` (orig says Char).
If a valid base is not provided, the function returns "N".

Examples  
≡≡≡≡≡≡≡≡≡≡

    julia> complement("A")
    "T"

    julia> complement("GCT")
    "CGA"

    julia> complement("AGCTACC")
    "TCGATGG"
"""
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

# ###5. Reverse Complement Function
"""
    reverse_complement(sequence)

    Takes a DNA sequence and returns the reverse complement
    of that sequence.
    
    Takes lowercase or uppercase sequences,
    but always returns uppercase.
    
    Examples
    ≡≡≡≡≡≡≡≡≡≡
        julia> reverse_complement("ACCTTT")
        "AAAGGT"
    
        julia> reverse_complement("CCGTAGTA")
        "TACTACGG"
    
        julia> rc = reverse_complement("ACAG");
    
        julia> println(rc)
        CTGT
    """
function reverse_complement(sequence)
    st= ""
    sequence =reverse(sequence)
    for i in 1:length(sequence)
        st = st*string(complement(sequence[i]))
    end
    return st
end

# ### Parse Function

"""
    function parse_fasta(path)

Reads a fasta-formated file and returns 2 vectors,
one containing parsed headers,
the other containing the entire sequence as a `String`.

Note: function does not validate DNA sequences for correctness.

Example
≡≡≡≡≡≡≡≡≡
    julia> ex1 = parse_fasta("data/ex1.fasta");

    julia> ex1 isa Tuple
    true

    julia> ex1[1]
    2-element Array{Tuple{String,String},1}:
      ("ex1.1", "easy")
      ("ex1.2", "multiline")

    julia> ex1[2]
    2-element Array{String,1}:
    "AATTATAGC"
    "CGCCCCCCAGTCGGATT"

    julia> ex2 = parse_fasta("data/ex2.fasta");

    julia> ex2[1]
    4-element Array{Tuple{String,String},1}:
      ("ex2.1", "oneper")
      ("ex2.2", "wrong")
      ("ex2.3", "wronger")
      ("ex2.4", "wrongest")
    
    julia> ex2[2]
    4-element Array{String,1}:
      "ATCCGT"
      "ATCGTGGaact"
      "ATCGTGGaact"
      "this isn't a dna string,but parse it anyway"
"""
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

# ### 7. Unique Kmer Function
"""
    function uniqueKmers(sequence, k)

Takes a sequence and a kmer length (k) 
and returns a list of strings of the unique kmers that appear within the DNA sequence.

Returns an "Invalid base" error if base is not within "AGCT".

Example
≡≡≡≡≡≡≡≡≡
    julia> uniqueKmers("ATGCGATG", 4)
    Set{Any} with 5 elements:
    "TGCG"
    "ATGC"
    "GATG"
    "CGAT"
    "GCGA"

    julia> uniqueKmers("ATGCN", 2)
    ERROR: Invalid base N encountered

    julia> uniqueKmers("ACT", 2)
    Set{Any} with 2 elements:
    "AC"
    "CT"
"""
function uniqueKmers(sequence, k)
    for base in sequence 
        if !occursin(base, "ACGT")
            error("Invalid base $base encountered")
        end
    end
    1 <= k <= length(sequence) || error("k must be a positive integer less than the length of the sequence")
    kmers = Dict()    
    stopindex = length(sequence) - k + 1
    ret = []
    for i in 1:stopindex
        kmer= sequence[i:i+k-1]
        kmer = uppercase(kmer)
        push!(ret, kmer)
        if haskey(kmers, kmer) == true
            kmers[kmer] = kmers[kmer] + 1
        else
            kmers[kmer] = 1
        end
    end
    return Set(ret)
end

# ### 8. KmerDistance Function
"""
    function kmerdist(set1, set2)

Takes two kmer sets and returns the distance between the kmer sets.
Must return a positive number between 0 and 1.

Returns 0.0 for identical sets.
Returns 1.0 for sets with no similarity.

Example
≡≡≡≡≡≡≡≡≡
    julia> kmerdist(uniqueKmers("ATCGATG",2), uniqueKmers("GCATACC",2))
    0.9

    julia> kmerdist(uniqueKmers("GCGCAT",2), uniqueKmers("ATAT",2))
    0.8

    julia> kmerdist(uniqueKmers("AATA", 4), uniqueKmers("CGGCCCG", 4))
    1.0
"""
function kmerdist(set1, set2)
    return 1 - (length(intersect(set1, set2))/length(union(set1,set2)))
end

# ###9. Kmertime Function
"""
    function kmertimedist(sequence1, sequence2)

Takes two kmer sets and returns the distance and time between the kmer sets.


Example
≡≡≡≡≡≡≡≡≡
    julia> kmertimedist(uniqueKmers("ATCGATG",2), uniqueKmers("GCATACC",2))
    0.9
    
    julia> kmertimedist(uniqueKmers("GCGCAT",2), uniqueKmers("ATAT",2))
    0.8

    julia> kmertimedist(uniqueKmers("AATA", 4), uniqueKmers("CGGCCCG", 4))
    1.0
"""

end

