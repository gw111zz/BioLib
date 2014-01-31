#!/usr/env/ruby -w

require 'test/unit'
require '../lib/biolib.rb'

class BioTest < Test::Unit::TestCase

  def test_read_fasta
      hash = BioLib.read_fasta('./fasta.txt')
      assert_equal({"Rosalind_0498"=>"AAATAAA", "Rosalind_2391"=>"AAATTTT", "Rosalind_2323" => "TTTTCCC", "Rosalind_0442"=>"AAATCCC", "Rosalind_5013"=>"GGGTGGGATATAAAATTTTTA"}, hash)
  end
 
  def test_count_nucleotides
    assert_equal({'A' => 0, 'C' => 0, 'G' => 0, 'T' => 0}, BioLib.tally_nucleotides(''))
    assert_equal({'A' => 3, 'C' => 0, 'G' => 0, 'T' => 0}, BioLib.tally_nucleotides('AAA'))
    assert_equal({'A' => 3, 'C' => 2, 'G' => 1, 'T' => 1}, BioLib.tally_nucleotides('CATACGA'))
    assert_equal({'A' => 20, 'C' => 12, 'G' => 17, 'T' => 21}, BioLib.tally_nucleotides('AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC'))
  end

  def test_transcribe_to_rna
    assert_equal('', BioLib.transcribe_to_rna('')) 
    assert_equal('GAUGGAACUUGACUACGUAAAUU', BioLib.transcribe_to_rna('GATGGAACTTGACTACGTAAATT')) 
  end

  def test_reverse_compliment
    assert_equal('', BioLib.reverse_compliment('')) 
    assert_equal('ACCGGGTTTT', BioLib.reverse_compliment('AAAACCCGGT'))
  end

  def test_gc_content
      assert_equal(0, BioLib.gc_content(''))
      assert_equal(0, BioLib.gc_content('ATTTTTTTTAAAAAAAATATATATTATATATAATATATATATATAATAT'))
      assert_equal(60.91954022988506, BioLib.gc_content('CCACCCTCGTGGTATGGCTAGGCATTCAGGAACCGGAGAACGCTTCAGACCAGCCCGGACTGGGAACCTGCGGGCAGTAGGTGGAAT'))
  end

  def test_hamming_distance
      same = 'ATGCATGCATGCATGCATGC'
      assert_equal(0, BioLib.hamming_distance(same, same))
      assert_equal(7, BioLib.hamming_distance('GAGCCTACTAACGGGAT', 'CATCGTAATGACGGCCT'))                   
  end

  def test_find_motif
      assert_equal([], BioLib.find_motif('', 'A'))
      assert_equal([], BioLib.find_motif('A', '' ))
      assert_equal([1], BioLib.find_motif('A', 'A' ))
      assert_equal([1], BioLib.find_motif('ATGC', 'AT' ))
      assert_equal([2], BioLib.find_motif('ATGC', 'TG' ))
      assert_equal([3], BioLib.find_motif('ATGC', 'GC' ))

      # Mutiple matches
      assert_equal([1,3,5,7,9], BioLib.find_motif('ATATATATAT', 'AT'))
      assert_equal([1,9], BioLib.find_motif('ATGCGCGCAT', 'AT'))
  end

  def test_find_adjacency_list
      sequences = { 's1' => 'AAATAAA', 's2' => 'AAATTTT', 's3' => 'TTTTCCC', 's4' => 'AAATCCC', 's5' => 'GGGTGGG' }
      assert_equal([ ['s1', 's2'], ['s1', 's4'], ['s2', 's3'] ], BioLib.adjacency_list(sequences, 1))
      assert_equal([ ['s1', 's2'], ['s1', 's4'], ['s2', 's3'] ], BioLib.adjacency_list(sequences, 2))
      assert_equal([ ['s1', 's2'], ['s1', 's4'], ['s2', 's3'] ], BioLib.adjacency_list(sequences, 3))

      sequences2 = { 's1' => 'AAAA', 's2' => 'AAAA' }
      assert_equal([], BioLib.adjacency_list(sequences2, 2))

      sequences3 = { 's1' => 'AAATAA', 's2' => 'AATAA' }
      assert_equal([ ['s1', 's2'], ['s2', 's1'] ], BioLib.adjacency_list(sequences3, 2))
  end

  def test_probability_phenotype
      assert_equal(0.7833333333333333, BioLib.probability_phenotype(2, 2, 2))
      assert_equal(0.7266949152542372, BioLib.probability_phenotype(15, 27, 18))
  end

  def test_num_permutations
      assert_equal(1, BioLib.num_permutations(1))
      assert_equal(120, BioLib.num_permutations(5))
  end

  def test_generate_permutations
      BioLib.generate_permutations(5)
  end

  def test_rna_to_protein
      assert_equal('MAMAPRTEINSTRING', BioLib.rna_to_protein('AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA'))
  end

  def test_compute_profile
      sequences = ['ATCCAGCT',
               	  'GGGCAACT',
             	    'ATGGATCT',
                  'AAGCAACC',
             	    'TTGGAACT',
             	    'ATGCCATT',
             	    'ATGGCACT']
     profile = BioLib.compute_profile(sequences)
     correct = {"A"=>[5, 1, 0, 0, 5, 5, 0, 0], "T"=>[1, 5, 0, 0, 0, 1, 1, 6], "G"=>[1, 1, 6, 3, 0, 1, 0, 0], "C"=>[0, 0, 1, 4, 2, 0, 6, 1]}
     assert_equal(correct, profile)

     assert_equal('ATGCAACT', (BioLib.compute_consensus(profile).join ''))
  end

  def test_calculate_protein_mass
      protein = 'SKADYEK'
      assert_equal(821.392, ('%.3f' % BioLib.calculate_protein_mass(protein)).to_f)
  end
end
