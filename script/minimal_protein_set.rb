#!/usr/bin/ruby

require 'yaml'
require 'set'
require 'optparse'
require 'ms/fasta'
require 'ms/fasta/ipi'

SET_RE = /Set\s+(.*)/i
PRECISION_EXT = ".precision.yml"

def sets_compare_to_paths(file, ext=PRECISION_EXT)
  dirname = File.dirname(File.expand_path(file))
  lines = IO.readlines(file).map {|v| v.chomp }.select {|v| v =~ /\w/}
  sets = {}
  current_set = nil
  lines.each do |line|
    if line =~ SET_RE
      current_set = $1.dup
      sets[current_set] = [] 
    else
      full_path = (File.join(dirname,(line + ext)))
      raise RuntimeError, "file #{full_path} does not exist!!" unless File.exist?(full_path)
      sets[current_set] << full_path
    end
  end
  sets
end

# returns those above the precision cutoff and those below
# requires a block for sorting the hits
def cut_by_precision(hits, mowse_index, prec_index, cutoff, &block)
  ###############################################
  ###############################################
  ###############################################
  ###############################################
  ###############################################
  # WORKING HERE
  ###############################################
  ###############################################
  ###############################################
  ###############################################
  low_to_hi_by_mowse = hits.sort_by &block

  found_lowest_above = false
  low_to_hi_by_mowse.partition do |hit|
    if ((!found_lowest_above) && (hit[prec_index] >= cutoff))
      found_lowest_above = true
    end
    found_lowest_above
  end
end

# returns [minimal_protein_to_uniq_peps_hash, indistinguishable_protein_hash] 
# takes a hash of proteins to aaseqs. Uses a greedy algorithm where
# things are sorted first by the number of uniq amino acid sequences and total
# aa length.  if a block is given, then will yield the prot and the
# peptide_array and sort by the returned value.  The greedy algorithm acts on
# the REVERSE of the sorted proteins.  indistinguishable_protein_hash is keyed
# on the proteins in the minimal_protein_array and gives an array of other
# proteins.  
def minimal_protein_set(proteins_to_aaseqs)
  blk_given = block_given?
  #STDERR.puts "using block for minimal_protein_set" if blk_given
  proteins_and_uniq_peps = []

  sorted_most_to_least = proteins_to_aaseqs.sort_by do |k,v| 
    if blk_given
      yield(k,v)
    else
      [ v.size, v.inject(0){|m,s| m+s.size} ]
    end
  end.reverse

  found_seq = Set.new

  same_peptide_hits = Hash.new{|h,k| h[k] = [] }

  last_peps = nil
  sorted_most_to_least.each do |prot, peps|
    sorted_peps = peps.sort # is it necessary to SORT?????????
    uniq_peps = peps.select do |pep|
      if found_seq.include?(pep)
        false
      else
        found_seq.add pep
        true
      end
    end
    if uniq_peps.size > 0
      proteins_and_uniq_peps << [prot, uniq_peps]
      last_peps = sorted_peps
    else
      if sorted_peps == last_peps 
        same_peptide_hits[proteins_and_uniq_peps.last.first] << prot
      end
    end
  end
  prot_to_uniq_peps_hash = {}
  proteins_and_uniq_peps.each do |prot, uniq_peps|
    prot_to_uniq_peps_hash[prot] = uniq_peps
  end

  [prot_to_uniq_peps_hash, same_peptide_hits]
end

def cutoffs_to_floats(ar, fdr=false)
  ar.map do |v|
    if v == 'nil' || v == '-'
      nil
    else
      answ = v.to_f
      fdr ? 1.0 - answ : answ
    end
  end
end

# returns a hash keyed on protein id that yields an array:
#   [#total hits, #aaseq+charge, #aaseq]
def stats_per_prot(prot_to_peps, seq_to_hits)
  per_protein_hash = {}
  prots_to_peps.each do |prot, uniq_pep_seqs|
    all_hits = Set.new
    aaseqcharges = Set.new
    aaseqs = Set.new

    uniq_pep_seqs.each do |pep_seq|
      all_hits = seq_to_hits[pep_seq]
      all.merge( all_hits )
      all_hits.each do |hit|
        aaseq = hit.aaseq
        aaseqs.add( aaseq )
        aaseqcharges.add( aaseq + hit.charge )
      end
      per_protein_hash[prot] = [all_hits.size, aaseqcharges.size, aaseqs.size]
    end
    per_protein_hash
end

opt = {
  :cutoffs => [nil],
  :outfile => "summary.yml",
}

opts = OptionParser.new do |op|
  op.banner = "usage: #{File.basename(__FILE__)} sets_compare.txt"
  op.separator "expects a sets_compare file in the format:"
  op.separator ""
  op.separator "Set <whatever>"
  op.separator "filename1_no_ext"
  op.separator "filename2_no_ext"
  op.separator "Set <whatever2>"
  op.separator "..."
  op.separator ""
  op.separator "appends to #{opt[:outfile]}, but will overwrite sections if the same precision_cutoff:"

  op.separator "- precision_cutoff: <float or nil>"
  op.separator "  num_unique_aaseqs: <int>"
  op.separator "  num_uniqe_aaseqs_charge: <int>"
  op.separator "  num_hits: <int>"
  op.separator ""
  op.separator "other information "
  op.separator "each file should already have a <file>.precision.yml file"
  op.separator ""
  op.on("-p", "--precision <0-1[,...]>", Array, "precision cutoff(s)") {|v| opt[:cutoffs] = cutoffs_to_floats(v)}
  op.on("-f", "--fdr <0-1[,...]>", Array, "false discovery rate cutoff (1-precision)") {|v| opt[:cutoffs] == cutoffs_to_floats(v, true) }
  op.separator ""
  op.separator "NOTE: if a precision cutoff is given, all hits that have a better"
  op.separator "score than the worst score at the cutoff are included, even if "
  op.separator "the precision for that hit was below the cutoff"
  op.separator "this prevents early, local aberrations in precision from messing"
  op.separator "up the analysis"
  op.separator ""
  op.on("--proteins <fasta>,<pep-db>", Array, "path to fasta and peptide centric DB", "peptide_centric_db is in the format: ", "<PEPTIDE>: <ID>-<ID>-<ID>") {|v| opt[:proteins] = v }
  op.on("--print-proteins", "gives an array of proteins") {|v| opt[:proteins] = v }
end

# later on we could implement full isoform resolution like IsoformResolver
# for now we will generate a report, realizing that some isoforms may not be
# reported
# it is implemented by using a pre-made map from sequence to protein groups
# then, a set of sequences allows one to deduce all the relationships from the 
# protein groups.

opts.parse!

if ARGV.size != 1
  puts opts.to_s
  exit
end


results = {}

protein_info = {}
results['protein_info'] = protein_info
results['results'] = []

sets_hash = sets_compare_to_paths(ARGV.shift)

if opt[:proteins]
  (fasta, pep_db_file) = opt[:proteins]
  
  # a hash indexed on ipi containing all info
  prot_header_hash = {}

  STDERR.print "Loading information from fasta file..."
  start = Time.now
  prot_sizes_hash = {}
  Ms::Fasta.open(fasta, 'rb', :io_index => []) do |obj|
    obj.each do |entry|
      hash = Ms::Fasta::Ipi.parse(entry.header)
      ipi = hash['IPI']
      prot_header_hash[ipi] = hash
      prot_sizes_hash[ipi] = entry.sequence.size
    end
  end
  STDERR.puts "#{Time.now - start} seconds."

  STDERR.print "Loading peptide centric DB (this takes about a minute)..."
  start = Time.now
  pep_db = YAML.load_file(pep_db_file)
  STDERR.puts "#{Time.now - start} seconds."


end

opt[:cutoffs].each do |cutoff|

  #results.reject {|hash| hash[:precision_cutoff] == cutoff } # clear out older results
  cutoff_results = {'precision_cutoff' => cutoff }
  results['results'] << cutoff_results

  #########################
  # FOR EACH SET:
  #########################
  sets_hash.each do |set, files|
    set_results = {}
    cutoff_results[set] = set_results

    # assumes the indices are the same into each data file

    # get the complete set of passing hits
    all_passing_hits = files.inject([]) do |all_passing_hits, file|
      hash = YAML.load_file(file)

      header_hash = hash['headers']
      Pep = Struct.new(*(header_hash.map {|v| v.to_sym }))
      passing_hits = 
        if cutoff
          (above, below) = cut_by_precision(hash['data'], mowse_index, prec_index, cutoff)
          above
        else
          hash['data']
        end
      all_passing_hits.push(*passing_hits)
    end

    
    # create an index from aaseq to hits
    seq_to_hits = Hash.new {|h,k| h[k] = []}
    uniq_seqcharge = Set.new
    all_passing_hits.each do |hit|
      seq_to_hits[hit[aaseq_index]] << hit
      uniq_seqcharge.add( hit[aaseq_index] + hit[charge_index] )
    end

    # determine the number of uniq aaseqs
    uniq_seqs = seq_to_hits.size

    ## THIS METHOD is probably faster, but much longer than just making the
    ## Set above.
    ##
    # determine the number of uniq aaseq+charge 
    #num_uniq_seqcharges = seq_to_hits.inject(0) do |sum, pair|
    #  (seq, hits) = pair
    #  set = Set.new
    #  initial = hits.first[charge_index]
    #  hits.each do |hit|
    #    charge = hit[charge_index]
    #    puts "FOUND ONE!!!!!!" if charge != initial
    #    set.add hit[charge_index]
    #    initial = charge
    #  end
    #  sum + set.size
    #end

    num_uniq_seqcharges = uniq_seqcharge.size

    set_results.merge!( { 'num peptide hits' => all_passing_hits.size,
      'num_uniq_aaseqs' => uniq_seqs,
      'num_uniq_aaseqs_charge' => num_uniq_seqcharges,
    })

    if opt[:proteins]

      # create an index from proteins to peptides
      prots_to_peps = Hash.new {|h,k| h[k] = [] }
      peptides_not_found = []
      seq_to_hits.keys.each do |seq|
        if pep_db.key?(seq)
          pep_db[seq].split('-').each do |prot|
            prots_to_peps[prot] << seq
          end
        else
          peptides_not_found << seq
        end
      end

      # Determine the number of 1) hits, 2) aaseqs, 3) aaseqcharges per protein BEFORE minimization


      # get the minimal protein set
      (prot_to_uniq_peps_hash, indistinguishable_protein_hash) = minimal_protein_set(prots_to_peps) do |prot,peps|
        # will sort with lowest 
        [ peps.size, peps.inject(0){|m,s| m+s.size}, -(prot_sizes_hash[prot])]
      end

      prot_to_uniq_peps_hash.each do |prot, peps|
        [prot, *indistinguishable_protein_hash[prot]].each do |prot|
          protein_info[prot] = prot_header_hash[prot]
        end
      end

      # some proteins in the file are reported as having an empty set of indistinguishable
      # proteins and I'm not sure why... delete them here..
      to_delete = []
      indistinguishable_protein_hash.each do |k,v|
        if v == []
          to_delete << k
        end
      end
      to_delete.each do |del| indistinguishable_protein_hash.delete(del) end

      set_results['indistinguishable_proteins'] = indistinguishable_protein_hash
      set_results['proteins'] = prot_to_uniq_peps_hash
      set_results['num proteins'] = prot_to_uniq_peps_hash.size
      set_results['num_aaseqs_not_in_pep_db'] = peptides_not_found.size
    end
  end
end

File.open(opt[:outfile], 'w') do |out|
  out.print results.to_yaml
end


