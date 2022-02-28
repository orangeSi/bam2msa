require "admiral"
require "./readfasta"

class Bam2Msa < Admiral::Command
  define_argument ref,
    description: "ref fasta file",
    required: true
  define_argument bam,
    description: "bam alignemnt file",
    required: true
  define_argument regions : String,
    description: "display read and ref msa alignment in these regions, example: chr1:1000-1200,chr2:2000-2300",
    required: true
  define_flag primary_only : Int32,
    default: 1_i32,
    description: "only for primary alignment. 0 mean all alingment, 1 is only primary alignment"
  define_flag span_whole_region_read_only : Int32,
    default: 1_i32,
    description: "only for read which span the whole region. 0 mean all read which overlap with the region, 1 mean is read which span the whole region"
  define_help description: "convert bam to msa format for alignment file"
  define_version "0.0.2"

  COMPILE_TIME = Time.local

  def run
    if ARGV.size == 0
      puts "Contact: ilikeorangeapple@gmail.com or go to https://github.com/orangeSi/bam2msa/issues"
      Bam2Msa.run "--help"
    end

    rgs = parser_regions(arguments.regions)
    #t0 = Time.utc
    ref = read_fasta(arguments.ref, chrs: rgs.keys)
    #t1 = Time.utc
    #puts "read_fasta cost #{t1-t0}s"
    #puts "flags.span_whole_region_read_only #{flags.span_whole_region_read_only}"
    #t2 = Time.utc
    arguments.regions.split(",").each do |rg|
      convert_bam2msa(arguments.bam, ref, flags.primary_only, rg, flags.span_whole_region_read_only)
    end
    #t3 = Time.utc
    #puts "convert_bam2msa cost #{t3-t2}s"
  end

  def convert_bam2msa(bam : String, ref : Hash(String, String), primary_only : Int32, region : String, span_whole_region_read_only : Int32)
    raise "error: cann't find samtools in $PATH" unless Process.find_executable("samtools")
    # title
    # puts "#refid\tref_cut_region\tref_msa\tquery_id\tquery_msa\tconsensus_msa\traw_cigar\tflag"
    puts "#query_msa\tref_msa\tconsensus_msa\tref_cut_region\tquery_id\tFLAG\tPOS_in_Bam\tCIGAR"

    # get region of output
    rgs = parser_regions(region)

    # read the bam
    # # check bam.bai file if exists!
    if File.exists?("#{bam}.bai") == false
      raise("error: #{bam}.bai not exists, so need samootls index #{bam} before!")
    end
    Process.run("samtools view #{bam} #{region}", shell: true) do |proc|
      while line = proc.output.gets
        # to do: parallel this by channel.send(line)
        msa = bam2msa_oneline(line, primary_only, ref, rgs, span_whole_region_read_only)
        next if msa.is_a?(Nil)
        # puts "#{msa.ref}\t#{msa.ref_region}\t#{msa.ref_msa}\t#{msa.query}\t#{msa.query_msa}\t#{msa.consensus}\t#{msa.cigar}\t#{msa.flag}"
        puts "#{msa.query_msa}\t#{msa.ref_msa}\t#{msa.consensus}\t#{msa.ref}:#{msa.ref_region}\t#{msa.query}\t#{msa.flag}\t#{msa.pos}\t#{msa.cigar}"
      end
    end
  end

  def bam2msa_oneline(line : String, primary_only : Int32, ref : Hash(String, String), rgs = {} of String => RG, span_whole_region_read_only : Int32 = 1)
    # puts "read line: #{line}"
    arr = line.split(/\t/)
    rflag = arr[1].to_i32
    ref_id = arr[2]

    return nil if !rgs.has_key?(ref_id)              # skip alignemtn which ref is not in regions
    return nil if primary_only && (rflag & 2304) > 0 # not primary alignment + supplementary alignment = 2304
    return nil if (rflag & 4) > 0                    # unmap = 4, filtered unmaped read

    query_id = arr[0]
    arr[0] = "_R_#{query_id}" if (rflag & 16) != 0 # reversed read map to ref

    # ref[arr[2]] # ref fasta
    # arr[9] # SEQ
    # arr[5] # cigar
    pos = arr[3].to_i32
    cigar = arr[5]
    raw_cigar = arr[5]

    rgs_start = rgs[ref_id].s
    rgs_end = rgs[ref_id].e
    the_region = "#{rgs_start}-#{rgs_end}"
    raise "error: ref #{ref_id} total size #{ref[ref_id].size}bp, but set #{rgs_end} in --regions \n" if rgs_end > ref[ref_id].size

    return nil if span_whole_region_read_only >=1 &&  pos > rgs_start

    # todo:if specific --regions chr1:200-300 , then cut the cigar for --regions
    # if rgs.has_key?(ref_id)
    # puts "raw_cigar=#{cigar} to #{ref_id}"
    # cigar, pos = cut_cigar_by_region(cigar, pos, rgs[ref_id].s, rgs[ref_id].e, query_id, ref_id)
    # return nil if cigar == ""
    # puts "new_cigar=#{cigar} to #{ref_id}" #sikaiwei
    # check the cigar form
    # raise "error:get wrong cigar=#{cigar} for #{query_id} #{ref_id}\n" if cigar =~ /[^\d^M^I^D^N^S^H^P^=^X]/
    # end

    ref_msa = ""
    raise "error: ref id #{ref_id} in the bam is not exists in the ref file\n" unless ref.has_key?(ref_id)
    ref_seq = ref[ref_id]
    ref_pos = -1 + pos

    query_msa = ""
    query_seq = arr[9]
    query_pos = 0

    # ---------- ref
    #    -------- query, the below code is for this
    # if ref_pos >= 1
    #  ref_msa = ref_seq[0...ref_pos]
    #  # (0...ref_pos).each { |e| query_msa += "-" }
    #  query_msa += "-"*ref_pos
    # end

    # trim the softclip in the start/end  of ref
    if cigar =~ /^(\d+)S/
      query_pos += $1.to_i32
      cigar = cigar.sub(/^\d+S/, "")
    end
    if cigar =~ /(\d+)S$/
      cigar = cigar.sub(/\d+S$/, "")
    end

    # loop the cigar
    cigar.scan(/\d+[MIDNSHP=X]/).each do |e|
      # puts "e=#{e}"
      if e[0] =~ /(\d+)(.)/
        clen = $1.to_i32
        ctype = $2 # MIDNSHP=X, reference to http://samtools.github.io/hts-specs/SAMv1.pdf
        # puts "#{clen}#{ctype}"
        case ctype
        when "M", "=", "X"
          ref_msa += ref_seq[ref_pos...ref_pos + clen]
          ref_pos += clen
          query_msa += query_seq[query_pos...query_pos + clen]
          query_pos += clen
          break if ref_pos > rgs_end # 跳过超出region的cigar
        when "I", "S"
          # (0...clen).each { |e| ref_msa += "-" } # add - to ref
          ref_msa += "-"*clen
          query_msa += query_seq[query_pos...query_pos + clen]
          query_pos += clen
        when "D", "N"
          ref_msa += ref_seq[ref_pos...ref_pos + clen]
          ref_pos += clen
          # (0...clen).each { |e| query_msa += "-" } # add - to query
          query_msa += "-"*clen
          break if ref_pos > rgs_end # 跳过超出region的cigar
        else
          raise("error: not recoginize cigar #{e}")
        end
      end
    end

    return nil if span_whole_region_read_only >=1 &&  (ref_pos+1) < rgs_end

    # --------- query
    #  ------------ ref
    if ref_seq.size < ref_pos
      raise("error: ref_seq.size #{ref_seq.size} <  ref_pos #{ref_pos}")
    end

    if query_msa.size != ref_msa.size
      raise("error: size not equal for #{arr[0]} and cigar=#{cigar}. #{query_msa.size} != #{ref_msa.size}")
    end

    # raise "error: ref_seq #{arr[2]} not equal, cigar=#{cigar}, ref_msa=#{ref_msa}\n" if ref_msa.gsub(/-/, "") != ref_seq

    #the_region = "1-#{ref[ref_id].size}"
    # # cut msa by --regions
    # puts "cut msa by --regions"
    consensus = ""
    ref_msa_cut = ""
    query_msa_cut = ""
    index = pos - 1

    # puts "query_msa=#{query_msa}"
    # puts "ref_msa=  #{ref_msa}"
    # puts "index=#{index} rgs_start=#{rgs_start} rgs_end=#{rgs_end}"

    # when pos > rgs_start
    if pos > rgs_start
      ref_msa_cut += ref_seq[rgs_start - 1...pos - 1]
      query_msa_cut += "-"*(pos - rgs_start)
    end

    ref_msa.each_char_with_index do |c, i|
      if index >= (rgs_start - 1) && index <= (rgs_end - 1)
        ref_msa_cut += c
        query_msa_cut += query_msa[i]
      end
      index += 1 if c != '-'
    end

    # when ref_pos < rgs_end
    if (ref_pos+1) < rgs_end
      ref_msa_cut += ref_seq[ref_pos...rgs_end]
      query_msa_cut += "-"*(rgs_end-ref_pos)
    end

    ref_msa = ref_msa_cut
    query_msa = query_msa_cut
    # get consensue sequence
    # puts "get consens"
    # puts "query_msa=#{query_msa}"
    # puts "ref_msa=  #{ref_msa}\n"
    query_msa.split(//).zip(ref_msa.split(//)) do |qe, re|
      next if qe == ""
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
    # puts "get consus done"

    if query_msa.size != consensus.size
      raise("error: consensus and query_msa size not equal for #{arr[0]} and #{cigar}. #{query_msa.size} != #{consensus.size}, consensus=#{consensus}, query_msa=#{query_msa}, ref_msa=#{ref_msa}")
    end
    if query_msa.size != ref_msa.size
      raise("error: ref_msa and query_msa size not equal for #{arr[0]} and #{cigar}. #{query_msa.size} != #{ref_msa.size}, query_msa=#{query_msa}, ref_msa=#{ref_msa}") 
    end
    # puts "read line done"
    return MSA.new(ref_id, the_region, ref_msa, arr[0], query_msa, consensus, raw_cigar, rflag, pos)
    # puts "#{arr[2]}\t#{ref_msa}\t#{arr[0]}\t#{query_msa}\t#{consensus}\t#{cigar}\t#{arr[1]}"
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
    # puts "cigars=#{cigar}, crs=#{crs}"
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
        # puts "#{cr.s}-#{cr.e} -> #{cr.v}#{cr.t}, rg_start=#{rg_start}, rg_end=#{rg_end}"
        if cr.t == "M" || cr.t == "=" || cr.t == "X" || cr.t == "D" || cr.t == "N"
          if rg_start >= cr.s && rg_start <= cr.e
            # puts "1"
            if rg_end > cr.e
              # puts "11"
              # # ----- 5M
              # #   -------- region
              new_cigar += "#{cr.v - (rg_start - cr.s)}#{cr.t}"
            else
              # puts "12"
              # # -------- 5M
              # #   ---- region
              new_cigar += "#{cr.v - (rg_start - cr.s) - (cr.e - rg_end)}#{cr.t}"
              break
            end
          elsif cr.s >= rg_start && cr.s <= rg_end
            if cr.e > rg_end
              # #         ----- 5M
              # #   -------- region
              new_cigar += "#{cr.v - (cr.e - rg_end)}#{cr.t}"
              break
            else
              # #         ----- 5M
              # #   ------------- region
              new_cigar += "#{cr.v}#{cr.t}"
            end
          elsif cr.e < rg_start
            # # ----- 5M
            # #            ------- region
            next
          elsif cr.s > rg_end
            # #                  ----- 5M
            # #  ------- region
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
    #raise "error: only support one region instead of #{regions}\n" if regions =~ /,/
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
    # puts "rgs=#{rgs}"
    return rgs
  end

  struct RG
    property s, e, v, t

    def initialize(@s : Int32, @e : Int32, @v : Int32 = 0, @t : String = "")
    end
  end

  struct MSA
    property ref, ref_region, ref_msa, query, query_msa, consensus, cigar, flag, pos

    def initialize(@ref : String, @ref_region : String, @ref_msa : String, @query : String, @query_msa : String, @consensus : String, @cigar : String, @flag : Int32, @pos : Int32)
    end
  end
end

Bam2Msa.run
