require "admiral"
require "./readfasta"

class Bam2Msa < Admiral::Command
  define_argument ref,
    description: "ref fasta file",
    required: true
  define_argument bam,
    description: "bam alignemnt file",
    required: true
  define_flag regions : String,
    default: "",
    description: "default display the msa of whole ref for every read. If not, set the specific regsion, ex: chr1:1000-3000"
  define_flag primary_only : Int32,
    default: 1_i32,
    description: "only for primary alignment. 0 mean all alingment, 1 is only primary alignment"
  define_flag display_left_softclip : Int32,
    default: 1_i32,
    description: "display softclip in the start of ref. 0 mean not display, 1 mean display"
  define_flag display_right_softclip : Int32,
    default: 1_i32,
    description: "display softclip in the end of ref. 0 mean not display, 1 mean display"
  define_help description: "convert bam to msa format for alignment file"
  define_version "0.0.1"

  COMPILE_TIME = Time.local

  def run
    if ARGV.size == 0
      puts "Contact: ilikeorangeapple@gmail.com or go to https://github.com/orangeSi/bam2msa/issues"
      Bam2Msa.run "--help"
    end

    rgs = parser_regions(flags.regions)
    ref = read_fasta(arguments.ref, chrs: rgs.keys)
    convert_bam2msa(arguments.bam, ref, flags.primary_only, flags.regions, display_left_softclip: flags.display_left_softclip, display_right_softclip: flags.display_right_softclip)
  end

  def convert_bam2msa(bam : String, ref : Hash(String, String), primary_only : Int32, regions : String = "", display_left_softclip : Int32 = 1, display_right_softclip : Int32 = 1)
    raise "error: cann't find samtools in $PATH" unless Process.find_executable("samtools")
    # title
    puts "#refid\tref_cut_region\tref_msa\tquery_id\tquery_msa\tconsensus_msa\traw_cigar\tflag"

    # get regions of output
    rgs = parser_regions(regions)
    rgs_key_size = rgs.keys.size

    # read the bam
    Process.run("samtools view #{bam}", shell: true) do |proc|
      while line = proc.output.gets
        # to do: parallel this by channel.send(line)
        msa = bam2msa_oneline(line,  primary_only, rgs_key_size, ref, rgs, display_left_softclip: display_left_softclip, display_right_softclip: display_right_softclip)
        next if msa.is_a?(Nil)
        puts "#{msa.ref}\t#{msa.ref_region}\t#{msa.ref_msa}\t#{msa.query}\t#{msa.query_msa}\t#{msa.consensus}\t#{msa.cigar}\t#{msa.flag}" 
      end
    end
  end

  def bam2msa_oneline(line : String, primary_only : Int32, rgs_key_size : Int32, ref : Hash(String, String), rgs = {} of String => RG, display_left_softclip : Int32 = 1, display_right_softclip : Int32 = 1)
    arr = line.split(/\t/)
    rflag = arr[1].to_i32
    ref_id = arr[2]
  
    return nil if rgs_key_size > 0 && !rgs.has_key?(ref_id) # skip alignemtn which ref is not in regions
    return nil if primary_only && (rflag & 2304) > 0 # not primary alignment + supplementary alignment = 2304
    return nil if (rflag & 4) > 0                    # unmap = 4, filtered unmaped read

    query_id = arr[0]
    arr[0] = "_R_#{query_id}" if (rflag & 16) != 0 # reversed read map to ref
    
    # ref[arr[2]] # ref fasta
    # arr[9] # SEQ
    # arr[5] # cigar
    pos = arr[3].to_i32
    cigar = arr[5]

    # todo:if specific --regions chr1:200-300 , then cut the cigar for --regions
    #if rgs.has_key?(ref_id)
      #puts "raw_cigar=#{cigar} to #{ref_id}"
      #cigar, pos = cut_cigar_by_region(cigar, pos, rgs[ref_id].s, rgs[ref_id].e, query_id, ref_id)
      #return nil if cigar == ""
      #puts "new_cigar=#{cigar} to #{ref_id}" #sikaiwei
      # check the cigar form
      #raise "error:get wrong cigar=#{cigar} for #{query_id} #{ref_id}\n" if cigar =~ /[^\d^M^I^D^N^S^H^P^=^X]/
    #end
    

    ref_msa = ""
    raise "error: ref id #{ref_id} in the bam is not exists in the ref file\n" unless ref.has_key?(ref_id)
    ref_seq = ref[ref_id]
    ref_pos = -1 + pos

    query_msa = ""
    query_seq = arr[9]
    query_pos = 0

    # ---------- ref
    #    -------- query, the below code is for this
    if ref_pos >= 1
      ref_msa = ref_seq[0...ref_pos]
      #(0...ref_pos).each { |e| query_msa += "-" }
      query_msa += "-"*ref_pos
    end

    # trim the softclip in the start/end  of ref
    if display_left_softclip <= 0 && cigar =~ /^(\d+)S/
        query_pos += $1.to_i32
        cigar = cigar.sub(/^\d+S/, "")
    end

    if display_right_softclip <= 0 && cigar =~ /(\d+)S$/
      cigar = cigar.sub(/\d+S$/, "") 
    end
   
    

    # loop the cigar
    cigar.scan(/\d+[MIDNSHP=X]/).each do |e|
      if e[0] =~ /(\d+)(.)/
        clen = $1.to_i32
        ctype = $2 # MIDNSHP=X, reference to http://samtools.github.io/hts-specs/SAMv1.pdf
        # puts "#{clen}#{ctype}"
        if ctype == "M" || ctype == "=" || ctype == "X"
          ref_msa += ref_seq[ref_pos...ref_pos + clen]
          ref_pos += clen
          query_msa += query_seq[query_pos...query_pos + clen]
          query_pos += clen
        elsif ctype == "I" || ctype == "S"
          #(0...clen).each { |e| ref_msa += "-" } # add - to ref
          ref_msa  += "-"*clen
          query_msa += query_seq[query_pos...query_pos + clen]
          query_pos += clen
        elsif ctype == "D" || ctype == "N"
          ref_msa += ref_seq[ref_pos...ref_pos + clen]
          ref_pos += clen
          #(0...clen).each { |e| query_msa += "-" } # add - to query
          query_msa +=  "-"*clen
        end
      end
    end

    # --------- query
    #  ------------ ref
    if ref_seq.size > ref_pos
      ref_msa += ref_seq[ref_pos..]
      #(ref_pos...ref_seq.size).each { |e| query_msa += "-" }
      query_msa += "-"*(ref_seq.size-ref_pos)
    elsif ref_seq.size < ref_pos
      raise("error: ref_seq.size #{ref_seq.size} <  ref_pos #{ref_pos}")
    end

    if query_msa.size != ref_msa.size
      raise("error: size not equal for #{arr[0]} and cigar=#{cigar}. #{query_msa.size} != #{ref_msa.size}")
    end

    raise "error: ref_seq #{arr[2]} not equal, cigar=#{cigar}, ref_msa=#{ref_msa}\n" if ref_msa.gsub(/-/, "") != ref_seq

    # get consensue sequence
    consensus = ""
    query_msa.split(//).zip(ref_msa.split(//)) do |qe, re|
      if qe == re
        if qe != "-" # match
          consensus += "="
        else
          raise "error: get both - for query=#{query_id} and ref=#{ref_id} for cigar=#{cigar}\n"
        end
      else
        if qe != "-" && re != "-" # mismatch
          consensus += "X"
        else           # indel
          if qe == "-" # deletion
            consensus += "D"
          else # insertion
            consensus += "I"
          end
        end
      end
    end

    if query_msa.size != consensus.size
      raise("error: consensus and query_msa size not equal for #{arr[0]} and #{cigar}. #{query_msa.size} != #{consensus.size}, consensus=#{consensus}")
    end

    #property ref, ref_msa, query, query_msa, consensus, cigar, flag

    the_region = "1-#{ref[ref_id].size}"
		## cut msa by --regions
	   if rgs.has_key?(ref_id)
          rgs_start = rgs[ref_id].s
          rgs_end = rgs[ref_id].e
          the_region = "#{rgs_start}-#{rgs_end}"
          ref_msa_cut = ""
          query_msa_cut = ""
          consensus_cut = ""
          index = 0

          raise "error: ref #{ref_id} total size #{ref[ref_id].size}bp, but set #{rgs_end} in --regions \n" if rgs_end > ref[ref_id].size

 				  ref_msa.each_char_with_index do |c,i|
            if index >= (rgs_start-1) && index <= (rgs_end-1)
              ref_msa_cut += c
              query_msa_cut += query_msa[i]
              consensus_cut += consensus[i]
            end
            index +=1 if c != '-'
          end
          ref_msa = ref_msa_cut
          query_msa = query_msa_cut
          consensus = consensus_cut
    end
    return MSA.new(ref_id, the_region, ref_msa, arr[0], query_msa, consensus, cigar, rflag)
    #puts "#{arr[2]}\t#{ref_msa}\t#{arr[0]}\t#{query_msa}\t#{consensus}\t#{cigar}\t#{arr[1]}"
  end
  
  def drive_into_cigar(cigar, pos) # caculate start-end of every type of in the cigar
    crs = [] of RG
    ref_pos = pos
    query_pos = 1
    cigar.scan(/\d+[MIDNSHP=X]/).each do |e|
      e[0] =~ /(\d+)(.)/
      clen = $1.to_i32
      ctype = $2 # MIDNSHP=X, reference to http://samtools.github.io/hts-specs/SAMv1.pdf
      if ctype == "M" || ctype == "=" || ctype == "X"
        crs << RG.new(ref_pos, ref_pos + clen - 1, clen, ctype)
        ref_pos += clen
      elsif ctype == "I" || ctype == "S"
        crs << RG.new(ref_pos - 1, ref_pos, clen, ctype)
      elsif ctype == "D" || ctype == "N"
        crs << RG.new(ref_pos, ref_pos + clen - 1, clen, ctype)
        ref_pos += clen
      end
    end
    #puts "cigars=#{cigar}, crs=#{crs}"
    return crs
  end

  def cut_cigar_by_region(cigar : String, pos : Int32, rg_start : Int32, rg_end : Int32, read_id : String, ref_id : String)
    new_cigar = ""
    new_pos = 1 # based 1, rg_start is also based 1
    
    # caculate start-end of every type of in the cigar, then check the below 5 situations.
    crs = drive_into_cigar(cigar, pos) # crs = [] of RG.new(start, end, type), type is MIDLSH=X
    
    # sitution6:
    # -------------------- ref
    #         ---- region
    #         ****
    # cigar ------------
    if rg_start >= crs[0].s && rg_end <= crs[-1].e
      crs.each do |cr|
        #puts "#{cr.s}-#{cr.e} -> #{cr.v}#{cr.t}, rg_start=#{rg_start}, rg_end=#{rg_end}"
        if cr.t == "M" || cr.t == "=" || cr.t == "X" || cr.t == "D" || cr.t == "N"
          if rg_start >= cr.s  && rg_start <= cr.e
            #puts "1"
            if rg_end > cr.e
            #puts "11"
  	          ## ----- 5M
    	        ##   -------- region
              new_cigar += "#{cr.v - (rg_start - cr.s) }#{cr.t}"
            else
            #puts "12"
              ## -------- 5M
              ##   ---- region
              new_cigar += "#{cr.v - (rg_start - cr.s) - (cr.e - rg_end) }#{cr.t}"
              break
            end
          elsif cr.s >= rg_start && cr.s <= rg_end
            if cr.e > rg_end
            ##         ----- 5M
            ##   -------- region
              new_cigar += "#{cr.v - (cr.e - rg_end) }#{cr.t}"
              break
            else
            ##         ----- 5M
            ##   ------------- region
              new_cigar += "#{cr.v}#{cr.t}"
            end

          elsif  cr.e < rg_start
           ## ----- 5M
           ##            ------- region
              next
          elsif cr.s > rg_end
            ##                  ----- 5M
            ##  ------- region
            next
          else
					  raise "error: this a bug, please issue to me~ read=#{read_id}, cigar=#{cigar}"
          end
        elsif cr.t == "I" || cr.t == "S"
          if rg_start <= cr.e
            new_cigar += "#{cr.v}#{cr.t}"
          end
        else
           raise("error: not recognize cigar #{cigar}")
        end 
      end
      new_cigar = new_cigar.sub(/I$/, "S")
      new_pos = rg_start
      return new_cigar, new_pos
    end

    # sitution4:
    # ref ----------------------
    # region ------
    #           ***
    #           ------ cigar
    if rg_start <= crs[0].s && rg_end >= crs[0].e
      return "", 1
    end


    # sitution1:
    # ref --------------------
    # region ----
    #        DDDD  ------- cigar
    if rg_end < pos
      new_pos = rg_start
      new_cigar = "#{rg_end - rg_start + 1}D"
      return new_cigar, new_pos
    end


    # sitution5:
    # ------------------------ ref
    #   -------------- region
    #       ********
    # cigar --------
    if crs[0].s >= rg_start && crs[-1].e <= rg_end
      return cigar, pos
    end

    # sitution2:
    #   -------------------------- ref
    #               ---- region
    # cigar ------  DDD, rg_start is also based 1D
    if crs[-1].e < rg_start
      new_pos = rg_start
      new_cigar = "#{rg_end - rg_start + 1}D"
      return new_cigar, new_pos
    end

    # sitution3:
    #   ------------------------- ref
    #          ------ region
    #          ***
    # cigar ------
    if rg_start >= crs[0].s && rg_start <= crs[-1].e
      return "", 1
    end

    # sitution4:
    # ref ----------------------
    # region ------
    #           ***
    #           ------ cigar
    if rg_start <= crs[0].s && rg_end >= crs[0].e
      return "", 1
    end


    # sitution4:
    # ref ----------------------
    # region ------
    #           ***
    #           ------ cigar
    if rg_start <= crs[0].s && rg_end >= crs[0].e
      return "", 1
    end

    raise "error: occur unexcept situation for read #{read_id} and ref #{ref_id}. when rg_start=#{rg_start},rg_end=#{rg_end} and position=#{pos},cigar=#{cigar}\n"
  end


  def parser_regions(regions : String)
    # chr1:100-300,chr3:500-800
    rgs = {} of String => RG
    return rgs if regions == ""
    raise "error: only support one region instead of #{regions}\n" if regions =~ /,/
    regions.split(/,/).each do |e|
      next if e == ""
      if e =~ /^([^:]+):(\d+)-(\d+)$/
        if $2.to_i <= $3.to_i
          rgs[$1] = RG.new($2.to_i, $3.to_i)
        else
          rgs[$1] = RG.new($3.to_i, $2.to_i)
        end
      else
        raise "error: not support region #{e}, example: chr1:100-300\n"
      end
    end
    return rgs
  end

  struct RG
    property s, e, v, t

    def initialize(@s : Int32, @e : Int32, @v : Int32 = 0, @t : String = "")
    end
  end

  struct MSA
    property ref, ref_region, ref_msa, query, query_msa, consensus, cigar, flag
    def initialize(@ref : String, @ref_region : String, @ref_msa : String, @query : String, @query_msa : String, @consensus : String, @cigar : String, @flag : Int32)
    end
  end
  
end

Bam2Msa.run
