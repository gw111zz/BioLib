#!/usr/env/ruby -w

module BioLib

    # Reads a file containing DNA sequences in FASTA format
    #
    # Suitable for short sequences as the entire file is read before this
    # method returns.
    #
    # Returns a Hash where the key is the sequence name and the value is the sequence.
    # 
    # Raises an ArgumentError if filename does not exist.
    #
    def self.read_fasta(filename)
        raise ArgumentError if !File.exist? filename
        result = Hash.new
        File.open(filename, 'r') do |file|
            name = ''
            dna = ''
            file.each_line do |line|
                if line =~ /^\s*>([^\n]+)$/
                    if dna != ''
                        result[name] = dna 
                        dna = ''
                    end
                    name = $1
                elsif line =~ /([ATGC]+)/
                    dna = dna + $1
                end
            end
            result[name] = dna
        end
        return result
    end

    # Returns a hash of the tally of symbols (characters) in the input string
    #
    # string::  The string to tally.
    #
    def self.tally(string)
        raise ArgumentError, 'Argument string should be a String' unless string.is_a? String

        symbol = lambda { |x| x }
        groups = string.chars.group_by(&symbol)
        groups.each { |k,v| groups[k] = v.length }
    end

    # Returns a hash of the tally of the nucleotides in an DNA string
    #
    # Ensures that the output always contains a count of 'A', 'C', 'G', 'T' bases
    #
    # string::  The string containing the nucleotide sequence. Ensure that the sequence does
    #           not consist of mixed cases. (The sequence should be "ACTGAC" and not "ACTGac".)
    #
    # Returns a Hash of the tallied nucleotides.
    #
    # For example,
    #
    #   Bio.tally_nucleotides('ACCG') 
    #
    # will return
    #
    #   { 'A' => 1, 'C' => 2, 'G' => 1, 'T' => 0 }
    #
    # Note that if there are non-base letters in the string, those will be tallied too.
    #
    # Raises ArgumentError if the argument string is not a String.
    #
    def self.tally_nucleotides(string)
        raise ArgumentError, 'Argument string should be a String' unless string.is_a? String

        tallied_string = tally(string)

        add_missing_symbol = lambda { |hash, symbol| hash[symbol] = 0 unless hash.has_key? symbol }

        ['A', 'C', 'G', 'T'].each do |symbol|
            add_missing_symbol.call(tallied_string, symbol)
        end 
        tallied_string
    end

    # 
    # Transcribes a DNA sequence into an RNA sequence. 
    #
    # string::  The DNA sequence that will be transcribed. The bases should be in capital letters.
    #
    # Returns the transcribed sequence as a String.
    #
    # Raises ArgumentError if the argument string is not a String. Invalid bases (ones not in 
    # the set "ACGT" are left intact in the string.
    #
    def self.transcribe_to_rna(string)
        raise ArgumentError, 'Argument string should be a String' unless string.is_a? String

        string.gsub 'T', 'U'
    end

    # 
    # Used to obtain the sequence of the complimentary strand of DNA given one strand.
    #
    # string::  The sequence to compute the complimentary for.
    #
    # Raises ArgumentError if the argument string is not a String. Invalid bases (ones not in 
    # the set "ACGT" are left intact in the string.
    #
    def self.reverse_compliment(string)
        raise ArgumentError, 'Argument string should be a String' unless string.is_a? String
        string.reverse.tr('ACGT', 'TGCA')
    end

    # 
    # Calculates the GC-content of a DNA sequence. Returns the value as a percentage.
    #
    # string::  The sequence to calculate the GC-content of.
    #
    # Raises ArgumentError if the argument string is not a String. 
    #
    def self.gc_content(string)
        raise ArgumentError, 'Argument string should be a String' unless string.is_a? String
        return 0 if string.length == 0

        tally = tally_nucleotides(string)
        (tally['G'] + tally['C']).to_f / string.length.to_f * 100.0
    end

    # Returns the Hamming distance between the two string s and t
    #
    # Raises ArgumentError if either s or t or both are not String or if
    # they are strings of different lengths.
    def self.hamming_distance(s, t)
        raise ArgumentError, 'Arguments s and t should both be a String' unless s.is_a?(String) && t.is_a?(String)
        raise ArgumentError, 'The lengths of s and t are not the same' if s.length != t.length
        
        # This is definitely not the best way to do this
        i = 0
        count = 0
        i.upto s.length do
            count = count + 1 if s[i] != t[i]
            i = i + 1
        end
        count
    end

    # Searches the given DNA sequence for the specified motif.
    #
    # Internally this uses the Rabin-Karp algorithm of substring finding.
    #
    # Returns the positions of the motif as an Array of integers where
    # 1 signifies the start of the sequence.
    #
    def self.find_motif(sequence, motif_to_find)
        raise ArgumentError, 'Arguments sequence and motif_to_find should both be Strings' unless sequence.is_a?(String) && motif_to_find.is_a?(String)
        
        result = []
        return result if sequence.length == 0 || motif_to_find.length == 0
        return result if motif_to_find.length > sequence.length
        return [1] if sequence.length == 1 && motif_to_find == sequence

        # Not too big otherwise computing a hash of the entire 1kb string will overflow
        # Fixnum
        prime = 8807
        hash = lambda do |string| 
            total = 0
            string.each_char do |n|
                total = total + n.ord * prime
            end
            total
        end
        rolling_hash = lambda { |previous_hash, new_value, previous_value| previous_hash + new_value.ord * prime - previous_value.ord * prime }
        
        # Calculate the hash of the motif
        hmotif = hash.call(motif_to_find)

        # Calculate the hash of the first n letters of the sequence where n = length of the motif
        hs = hash.call(sequence[0, motif_to_find.length])

        i = 0
        while i < (sequence.length - motif_to_find.length) do
            if hmotif == hs
                result << (i + 1) if motif_to_find == sequence[i, motif_to_find.length]
            end
            hs = rolling_hash.call(hs, sequence[i + motif_to_find.length], sequence[i])
            i = i + 1
        end

        # Would be nice to not need this bit
        if hmotif == hs 
            result << (i + 1) if motif_to_find == sequence[i, motif_to_find.length]
        end

        result
    end

    # Computes the On overlap graph for the specified DNA sequences.
    #
    # sequences :: A Hash of sequences where the keys are the sequence name and the values
    #              are the sequences.
    #
    # o_value   :: The order of the graph to calculate.
    #
    # Returns an Array of Arrays of pairs of sequence names. For example, if the input sequences list
    # was { 's1' => 'AATT', 's2' => 'TTCC', 's3' => 'CCGG' }, to indicate that there are edges from
    # s1 -> s2 and s2 -> s3, the output would be [ ['s1', 's2'], ['s2', 's3'] ] 
    #
    def self.adjacency_list(sequences, o_value)
        raise ArgumentError, 'Argument sequence should be a Hash' unless sequences.is_a?(Hash)
        raise ArgumentError, 'Argument o_value should be a positive integer' unless o_value.is_a?(Fixnum) && o_value > 0

        is_directed_edge = lambda { |s, t, o| s != t && t.start_with?(s[-o, o]) }

        result = []
        names = sequences.keys

        for i in 0..names.length - 1
            for j in 0..names.length - 1
                 next if i == j
                 if is_directed_edge.call(sequences[names[i]], sequences[names[j]], o_value)
                    result << [names[i], names[j]]
                 end
            end
        end
        result
    end
    
    # Calculate the probability that the dominant phenotype is displayed in a population
    # of k + m + n organisms assuming any can mate with any other.
    #
    # k :: Number of homozygous organisms in the population
    # m :: Number of heterozygous organisms in the population
    # n :: Number of homozygous recessive organisms in the population
    #
    # Returns the probability that the dominant phenotype is displayed.
    #
    def self.probability_phenotype(k, m, n)
        raise ArgumentError, "k should be a positive integer" unless k > 0
        raise ArgumentError, "m should be a positive integer" unless m > 0
        raise ArgumentError, "n should be a positive integer" unless n > 0

        total = k + m + n
        kk = k.to_f / total.to_f * (k.to_f - 1.0) / (total - 1).to_f
        km = k.to_f / total.to_f * m.to_f / (total - 1).to_f
        mk = m.to_f /  total.to_f * k.to_f / (total - 1).to_f
        kn = k.to_f / total.to_f * n.to_f / (total - 1).to_f
        nk = n.to_f / total.to_f * k.to_f / (total - 1).to_f
        mm = (m.to_f / total.to_f * (m.to_f - 1.0) / (total - 1).to_f) * 0.75
        mn = (m.to_f / total.to_f * n.to_f / (total - 1).to_f) * 0.5
        nm = (n.to_f / total.to_f * m.to_f / (total - 1).to_f) * 0.5
        
        kk + km + mk + kn + nk + mm + mn + nm
    end

    # Returns the number of permutations of n integers.
    #
    # n:: A positive integer.
    #
    # For example, is n = 5, it will return all permutations of
    # 1,2,3,4,5
    def self.num_permutations(n)
        raise ArgumentError, "n should be a positive integer" unless n > 0
        return 1 if n == 1
        n.downto(1).reduce(:*)
    end

    # Returns an Array of all permutations of n in lexicographical
    # order.
    #
    # n:: The permutations of n will be returned. n should be a positive integer.
    #
    def self.generate_permutations(n)
       raise ArgumentError, "n should be a positive integer" unless n > 0
       return [1] if n == 1

       result = []

       permutation = []
       1.upto(n) { |i| permutation << i } 

       result << permutation

       next_permutation = lambda do |current|
           copy = current.dup
           i = copy.length - 2

            until copy[i] < copy[i+1]
                i = i - 1
            end

            return [] if i < 0

            j = copy.length - 1
            until copy[j] > copy[i]
                j = j - 1
            end
            copy[i], copy[j] = copy[j], copy[i]
            copy = copy[0, i + 1] + copy[i + 1, copy.length - i].reverse
            return copy
       end

       begin 
          new_perm = next_permutation.call(result[result.length - 1])
          result << new_perm unless new_perm == []
       end until new_perm == []

       result
    end

    # A table indicating the translation of individual RNA codons into amino acids for the purpose of protein creation.
    # '' means a stop
    CodonTable = { 'UUU' => 'F',   'CUU' => 'L', 'AUU' => 'I', 'GUU' => 'V',
                   'UUC' => 'F',   'CUC' => 'L', 'AUC' => 'I', 'GUC' => 'V',
                   'UUA' => 'L',   'CUA' => 'L', 'AUA' => 'I', 'GUA' => 'V',
                   'UUG' => 'L',   'CUG' => 'L', 'AUG' => 'M', 'GUG' => 'V',
                   'UCU' => 'S',   'CCU' => 'P', 'ACU' => 'T', 'GCU' => 'A',
                   'UCC' => 'S',   'CCC' => 'P', 'ACC' => 'T', 'GCC' => 'A',
                   'UCA' => 'S',   'CCA' => 'P', 'ACA' => 'T', 'GCA' => 'A',
                   'UCG' => 'S',   'CCG' => 'P', 'ACG' => 'T', 'GCG' => 'A',
                   'UAU' => 'Y',   'CAU' => 'H', 'AAU' => 'N', 'GAU' => 'D',
                   'UAC' => 'Y',   'CAC' => 'H', 'AAC' => 'N', 'GAC' => 'D',
                   'UAA' => '',    'CAA' => 'Q', 'AAA' => 'K', 'GAA' => 'E',
                   'UAG' => '',    'CAG' => 'Q', 'AAG' => 'K', 'GAG' => 'E',
                   'UGU' => 'C',   'CGU' => 'R', 'AGU' => 'S', 'GGU' => 'G',
                   'UGC' => 'C',   'CGC' => 'R', 'AGC' => 'S', 'GGC' => 'G',
                   'UGA' => '',    'CGA' => 'R', 'AGA' => 'R', 'GGA' => 'G',
                   'UGG' => 'W',   'CGG' => 'R', 'AGG' => 'R', 'GGG' => 'G'  }

    # Transcripts the supplied RNA string to a protein string
    #
    # rna:: The RNA string
    #
    def self.rna_to_protein(rna)
        raise ArgumentError, 'Argument sequence should be a string' unless rna.is_a?(String)
        raise ArgumentError, 'Argument sequence length should be a multiple of 3' unless rna.length % 3 == 0

        protein = ''
        (0..rna.length).step(3) do |chunk|
            key = rna[chunk, 3]
            raise ArgumentError, "Unknown rna sequence found #{key}" unless CodonTable.has_key? key
            return protein if CodonTable[key] == ''
            protein = protein + CodonTable[key] 
        end

        protein
    end

    # Computes the profile of the DNA sequences.
    #
    # sequences:: An Array of DNA sequences.
    #
    # Returns a Hash containing the tally of each base.
    def self.compute_profile(sequences)
        raise ArgumentError, 'Argument sequences should be an Array' unless sequences.is_a?(Array)

        tally = Hash.new
        tally['A'] = []
        tally['T'] = []
        tally['G'] = []
        tally['C'] = []

        sequences.each do |sequence|
            0.upto(sequences[0].length - 1) do |i|
                tally['A'] << 0 if !tally.has_key? 'A' or tally['A'].length <= i
                tally['T'] << 0 if !tally.has_key? 'T' or tally['T'].length <= i
                tally['G'] << 0 if !tally.has_key? 'G' or tally['G'].length <= i
                tally['C'] << 0 if !tally.has_key? 'C' or tally['C'].length <= i
                char = sequence[i]
                tally[char][i] = tally[char][i] + 1
            end
        end
        tally
    end

    # Computes the consensus of the DNA profile.
    #
    # profile:: A Hash with a base as the key and the tally as the value.
    #           A value for each of A, T, G and C should be present in the Hash.
    #
    # Returns an Array containing the consensus.
    #
    def self.compute_consensus(profile)
        raise ArgumentError, 'Argument profile should be a Hash' unless profile.is_a?(Hash)
        raise ArgumentError, 'Argument profile should have key A' unless profile.has_key?('A')
        raise ArgumentError, 'Argument profile should have key T' unless profile.has_key?('T')
        raise ArgumentError, 'Argument profile should have key G' unless profile.has_key?('G')
        raise ArgumentError, 'Argument profile should have key C' unless profile.has_key?('C')

        # todo check profile has right length
        consensus = []
        0.upto(profile['A'].length - 1) do |i|
            #puts "i = #{i}, T = #{profile['T'][i]}, G = #{profile['G'][i]}, C = #{profile['C'][i]}, A = #{profile['A'][i]}"
            consensus[i] = 'A'
            consensus[i] = 'T' if profile['T'][i] > profile[consensus[i]][i]
            consensus[i] = 'G' if profile['G'][i] > profile[consensus[i]][i]
            consensus[i] = 'C' if profile['C'][i] > profile[consensus[i]][i]
        end

        consensus
    end

    # The monoisotopic mass table for amino acids is a table listing the mass of each possible amino acid residue, 
    # where the value of the mass used is the monoisotopic mass of each residue.
    MonoisotopicMassTable = 
        { 'A' => 71.03711,  'C' => 103.00919, 'D' => 115.02694, 'E' => 129.04259, 'F' => 147.06841, 'G' => 57.02146,
          'H' => 137.05891, 'I' => 113.08406, 'K' => 128.09496, 'L' => 113.08406, 'M' => 131.04049, 'N' => 114.04293,
          'P' => 97.05276,  'Q' => 128.05858, 'R' => 156.10111, 'S' => 87.03203,  'T' => 101.04768, 'V' => 99.06841,
          'W' => 186.07931, 'Y' => 163.06333 }

    # Calculates the mass of the supplied protein.
    #
    # Returns the weighted string.
    def self.calculate_protein_mass(protein)
        raise ArgumentError, "protein should be a string" unless protein.is_a?(String) 

        mass = 0;
        protein.each_char do |c|
          raise ArgumentError, "Invalid char in protein string" if !MonoisotopicMassTable.has_key?(c)
          mass += MonoisotopicMassTable[c]
        end

        mass
    end
end


