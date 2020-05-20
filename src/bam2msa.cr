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
    description: "default display the msa of whole ref for every read. If not, set the specific regsion, ex: chr1:100-300,chr3:500-800"
  define_flag primary_only : Int32,
    default: 1_i32,
    description: "only for primary alignment. 0 mean all alingment, 1 is only primary alignment"
  define_flag display_softclip : Int32,
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
    msa = convert_bam2msa(arguments.bam, ref, flags.primary_only, flags.regions)
  end
  def convert_bam2msa(bam : String, ref : Hash(String, String), primary_only : Int32, regions : String = "")
    raise "error: cann't find samtools in $PATH" unless Process.find_executable("samtools")
    # title
    puts "#refid\tref_msa\tquery_id\tquery_msa\tconsensus_msa\tcigar\tflag"

    # get regions
    rgs = parser_regions(regions)
    rgs_key_size = rgs.keys.size

    Process.run("samtools view #{bam}", shell: true) do |proc|
      while line = proc.output.gets
        #puts line
        # qid, flag, refid, 
        arr = line.split(/\t/)
        rflag = arr[1].to_i32
        ref_id = arr[2]
        next if rgs_key_size > 0 && ! rgs.has_key?(ref_id)
        next if primary_only && (rflag & 256) != 0
        query_id = arr[0]
        arr[0] = "_R_#{query_id}" if (rflag & 16) != 0 # reversed read map to ref
        #ref[arr[2]] # ref fasta
        #arr[9] # SEQ
        #arr[5] # cigar
        pos = arr[3].to_i32
        cigar = arr[5]

        # if specific --regions chr1:200-300 , then cut the cigar for --regions
        if rgs.has_key?(ref_id)
          pos, cigar = cut_cigar_by_region(cigar, pos, rgs[ref_id].s, rgs[ref_id].e, query_id, ref_id)
        end
        
        ref_msa = ""
        ref_seq = ref[ref_id]
        ref_pos = -1 + pos

        query_msa = ""
        query_seq = arr[9]
        query_pos = 0
        
        if ref_pos != 0
          ref_msa = ref_seq[0...ref_pos] 
          (0...ref_pos).each {|e| query_msa += "-"}
        end
        #puts "ref_seq is #{ref_seq}, query_seq is #{query_seq}"  

        split_cigar(cigar).each do |e|
          if e =~ /(\d+)(.)/
            clen = $1.to_i32
            ctype = $2  # MIDNSHP=X, reference to http://samtools.github.io/hts-specs/SAMv1.pdf
            #puts "#{clen}#{ctype}"  
            if ctype == "M" || ctype == "=" || ctype == "X"
              ref_msa += ref_seq[ref_pos...ref_pos+clen]
              ref_pos += clen
              query_msa += query_seq[query_pos...query_pos+clen]
              query_pos += clen
            elsif ctype == "I" || ctype == "S" 
              (0...clen).each {|e| ref_msa +="-" } # add - to ref
              query_msa += query_seq[query_pos...query_pos+clen]
              query_pos += clen
            elsif ctype == "D" || ctype == "N" 
              ref_msa += ref_seq[ref_pos...ref_pos+clen]
              ref_pos += clen
              (0...clen).each {|e| query_msa +="-" } # add - to query
            end
          end
        end
        #puts ">#{arr[2]}\n#{ref_msa}\n>#{arr[0]}\n#{query_msa}"
        if ref_seq.size >= ref_pos
          ref_msa += ref_seq[ref_pos..]
          (ref_pos...ref_seq.size).each {|e| query_msa +="-"}
        else
          raise("error: ref_seq.size #{ref_seq.size} ref_pos #{ref_pos}")
        end
        if query_msa.size != ref_msa.size
          raise("error: size not equal for #{arr[0]} and #{cigar}. #{query_msa.size} != #{ref_msa.size}")
        end
        
        raise "error: ref_seq #{arr[2]} not equal\n" if ref_msa.gsub(/-/, "") != ref_seq

        # get consensue sequence
        consensus = ""
        query_msa.split(//).zip(ref_msa.split(//)) do |qe, re|
          if qe == re
            if qe != "-" # match
              consensus +="="
            else
              raise "error: get both - for query=#{query_id} and ref=#{ref_id} for cigar=#{cigar}\n"
            end
          else
            if qe != "-" && re != "-" #mismatch
              consensus +="X"
            else # indel
              if qe == "-" # deletion
                consensus +="D"
              else # insertion
                consensus +="I"
              end
              
            end
          end
        end
          
        if query_msa.size != consensus.size
          raise("error: consensus and query_msa size not equal for #{arr[0]} and #{cigar}. #{query_msa.size} != #{consensus.size}, consensus=#{consensus}")
        end
        
        puts "#{arr[2]}\t#{ref_msa}\t#{arr[0]}\t#{query_msa}\t#{consensus}\t#{cigar}\t#{arr[1]}"
        
        
      end
    end
  end
  
  def drive_into_cigar(cigar, pos) #caculate start-end of every type of in the cigar
    
    
    
    return crs
  end

  def cut_cigar_by_region(cigar : String, pos : Int32, rg_start : Int32, rg_end : Int32, read_id : String, ref_id : String)
    new_cigar = ""
    new_pos = 0
    
    #sitution1:
      #region ----  
      #       DDDD  ------- cigar
    if rg_end < pos
      new_pos = rg_start
      new_cigar = "#{rg_end - rg_start +1}D"
      return new_cigar, new_pos
    end

    
   # caculate start-end of every type of in the cigar, then check the below 5 situations.
    crs = drive_into_cigar(cigar, pos) # crs = [] of RG.new(start, end, type), type is MIDLSH=X
    
    #sitution2:
      #              ---- region
      #cigar ------  DDDD
    if crs[-1].e < rg_start
      new_pos = rg_start
      new_cigar = "#{rg_end - rg_start +1}D"
      return new_cigar, new_pos
    end

    #sitution3:
      #          ------ region
      #          **
      #cigar ------ 
    
    #sitution4:
      # region ------ 
      #           ***
      #           ------ cigar
    
    #sitution5:
      #   -------------- region
      #      ********
      #cigar --------
    
    #sitution6:
      #         ---- region
      #         ****
      #cigar ------------
    
    raise "error: occur unexcept situation for read #{read_id} and ref #{ref_id}. when rg_start=#{rg_start},rg_end=#{rg_end} and position=#{pos},cigar=#{cigar}\n"  
    
  end

  def parser_regions(regions : String)
    #chr1:100-300,chr3:500-800
    rgs = {} of String => RG
    return rgs if regions == ""
    regions.split(/,/).each do |e|
      next if e == ""
      if e =~ /^([^:]+):(\d+)-(\d+)$/
        if $2 <= $3
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
    property s, e, v
    def initialize(@s : Int32, @e : Int32, @v : String = "")
    end
  end

  def split_cigar(cigar : String)
  clen = ""
  ctype = ""
  cigars = [] of String
  cigar.split(//).each do |e|
    if e =~ /[MIDNSHP=X]/
      ctype = e
      clen = clen.to_i32
      cigars << "#{clen}#{ctype}"
      clen = ""
      ctype = ""
    elsif e=~ /[0-9]/
      clen = "#{clen}#{e}"
    else
      raise("error: get wrong cigar #{e} in #{cigar}")
    end
  end
  return cigars
  end


end

Bam2Msa.run
