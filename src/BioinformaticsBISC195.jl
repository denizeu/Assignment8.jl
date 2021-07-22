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
export kmerloc
export kmertime
export kmertimes

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
    
    julia> gc_content("GGCCNNNN")
    1.0
"""
function gc_content(sequence)
    comp = composition(sequence)
    ng = get(comp, 'G', 0) #etc
    nc = get(comp, 'C', 0)
    na= get(comp, 'A', 0)
    nt= get(comp, 'T', 0)
    return (ng + nc) / (ng + nc + nt + na) #divided by sum of AGCT w/ composition function
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

# ### Mean Lengths and Counts
"""
    lengthcount(path)

    Takes data and returns the mean and standard deviation of sequence lengths and of gc content.
    
    Takes lowercase or uppercase sequences,
    but always returns uppercase.
    
    Examples
    ≡≡≡≡≡≡≡≡≡≡
        julia> lengthcount("data/example.fasta")
        (44806.2566948847, 45.93354112164542, 0.7767768512386324, 0.00047821765440753630)

        julia> lengthcount("data/example.fasta")
        (25714.4156356907, 68.78024813244301, 0.3631206332359013, 0.00021470921328415272)
    """
function lengthcount(path)
    lengths= []
    counts= []
    data= parse_fasta(path)
    for i in data[2]
        push!(lengths, length(i))
        gs = count(==('G'), i)
        cs = count(==('C'), i)
        GCcontent= (gs + cs)/length(i)
        push!(counts, GCcontent)
    end
    return (mean(lengths), std(lengths), mean(counts), std(counts))
end

# ### Minimum and Maximum Function
"""
    minMax(path)

    Takes data and returns the minimum and maximum of sequence lengths.
    
    Examples
    ≡≡≡≡≡≡≡≡≡≡
        julia> minMax("data/ex1.fasta")
        (9, 17)
    """
function minMax(path)
    ret = parse_fasta(path)[2]
    retmin = length(minimum(ret))
    retmax = length(maximum(ret))
    return (retmin, retmax)
end

# ### Sequence Lengths Function
"""
    lengthcount(path)

    Takes data and returns the sequence lengths.

    Examples
    ≡≡≡≡≡≡≡≡≡≡
        julia> seqlength("data/ex1.fasta")
        9
        17
    """
function seqlength(path)
    lengths = []
    data = parse_fasta(path)
    for i in data[2]
         push!(lengths, length(i))
    end
    return lengths
end

# ### Unique Kmer Function
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
        if !occursin(base, "ACGTN")
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
    return collect(Set(ret))
end
# ### Sorting Sequences to 5%
"""
    function sorting(path)

Takes a path and returns a vector with the sequence lengths greater than 29500.

Example
≡≡≡≡≡≡≡≡≡
    julia> sorting(path)
    Set{Any} with 5 elements:

"""
function sorting(path)
    sequences = seqlength(path)
    headers = parse_fasta(path)[1]
    check = findall(x->x<29500, sequences)
    deleteat!(sequences, check)
    deleteat!(headers, check)
    return (sequences, headers)
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
end


 # ###9. Kmertime Function
"""
    function kmertime(path)

Takes a dataset and returns three arrays of the unique kmers for early, middle, and late time periods.

Example
≡≡≡≡≡≡≡≡≡
    julia> kmertime("data/fake_sample_dataset.fasta")
    early=["ACAT", "AGAT", "CCAT"]
    middle= ["GGTA", "AAAC", "ATAC"]
    late= ["CGAT", "ACCA", "TATC"]
"""
function kmertime(path) #whatto set k to?
    k = 40
    early= []
    middle= []
    late= [] #Initializing empty arrays for 3 time periods: early(2019), middle(2020), and late(2021)
    headers = parse_fasta(path)[1]  
    sequences = parse_fasta(path)[2] 
    for i in 1:length(sequences)
        if occursin("2019", headers[i]) #Findall occurrences of 2019 per header
            push!(early, uniqueKmers(sequences[i], k)) #If the header contains the date "2019", the kmer will be pushed into the "early" array, calling the uniqueKmer function to process how many unique kmers exist in the sequence
        end
        if occursin("2020", headers[i])
            push!(middle, uniqueKmers(sequences[i], k)) #If the header contains the date "2020", the kmer will be pushed into the "middle" array
        end
        if occursin("2021", headers[i])
            push!(late, uniqueKmers(sequences[i], k)) #If the header contains the date "2021", the kmer will be pushed into the "late" array
        end
    end
    return (early, middle, late)
end

function kmertimes(path)
    early = length(kmertime(path)[1])
    middle = length(kmertime(path)[2])
    late = length(kmertime(path)[3])
    cptearly = 0 
    cptmiddle = 0 
    cptlate = 0
    for i in 1:early 
        cptearly = cptearly + 1
    end
    for i in 1:middle
        cptmiddle  = cptmiddle + 1
    end
    for i in 1:late
        cptlate  = cptlate + 1
    end
    return(cptearly, cptmiddle, cptlate)
end

### Kmertime Plot
using Plots #x=time period (early, middle late), y= number of unique kmers of size "k"
histogram(kmertimes("data/genomes_CoV2.fasta")) #creating a histogram for the time periods vs. number of unique kmers

# ###10. Kmerlocation Function
"""
    function kmerloc(path)

Takes a dataset and returns two arrays of the unique kmers for Turkey and Japan.
Then, it finds the distance between the two largest kmer sets in Turkey and Japan.

Example
≡≡≡≡≡≡≡≡≡
    julia> kmerloc("data/fake_sample_dataset.fasta")
    Tu = "ACAT", "GCGT", "AGTA"
    Japan= "CCGT", "ATAT", "GGGC", "ACGC"
    0.6 
"""
function kmerloc(path)
    Tu= []
    Ja= []
    for header in eachline(path)
        if header occursin("TUR", path) #if TUR is located in the header, push the unique kmers to Turkey's array
            push!(Tu, uniqueKmers(sequence, k))
        end
        if header occursin("Japan", path)#if Japan is located in the header, push the unique kmers to Japan's array
            push!(Ja, uniqueKmers(sequence, k))
        end
            return Tu, Ja #return the two arrays
        kmerdist(Tu, Ja) #find the distance between the greatest kmerset in each (--How do i locate the greatest of each??)
        end
        return kmerdist #returns the distance between the two for the two countries
    end

 ### Kmer Location Plot
Plots.gr()
x= ["Turkey", "Japan"] #x-value is location: turkey or japan
y= [Tu, Ja] #y value holds the arrays of unique kmers 
pie(x, y, title= "Number of Unique COVID-19 Kmers in Turkey vs. Japan") #piechart shows highest number of unique kmers for each location=#

