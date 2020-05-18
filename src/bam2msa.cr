require "admiral"
class Bam2Msa < Admiral::Command
	define_argument ref,
		description: "ref fasta file",
		required: true
	define_argument bam,
		description: "bam alignemnt file",
		required: true
	define_flag primary_only : Int32,
		default: 1_i32,
		description: "only for primary alignment. 0 mean all alingment, 1 is only primary alignment"
	define_help description: "convert bam to msa format"
	define_version "0.0.1"

	COMPILE_TIME = Time.local

	def run
		if ARGV.size == 0
			puts "Contact: ilikeorangeapple@gmail.com or go to https://github.com/orangeSi/grepfile/issues"
			Bam2Msa.run "--help"
		end
		ref = read_fasta(arguments.ref)
		msa = convert_bam2msa(arguments.bam, ref, flags.primary_only)
	end
	def read_fasta(f : String)
		ref = {} of String => String
		return  ref
	end
	def convert_bam2msa(bam : String, ref : Hash(String, String), primary_only : Int32)
		puts "#refid\tref_msa\tquery_id\tquery_msa\tcigar\tflag"
		Process.run("samtools view #{bam}", shell: true) do |proc|
			while line = proc.output.gets
				#puts line
				# qid, flag, refid, 
				arr = line.split(/\t/)
				next if primary_only && (arr[1].to_i32 & 256) != 0 
				arr[0] = "_R_#{arr[0]}" if (arr[1].to_i32 & 16) != 0 # reversed read map to ref
				#ref[arr[2]] # ref fasta
				#arr[9] # SEQ
				#arr[5] # cigar
				ref_msa = ""
				ref_seq = ref[arr[2]]
				query_msa = ""
				query_seq = arr[9]
				ref_pos = 0
				query_pos = 0
				cigar = arr[5]
				#puts "ref_seq is #{ref_seq}, query_seq is #{query_seq}"	

				split_cigar(cigar).each do |e|
					if e =~ /(\d+)(.)/
						clen = $1.to_i32
						ctype = $2  # MIDNSHP=X
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
				if query_msa.size != ref_msa.size
					raise("error: size not equal for #{arr[0]} and #{cigar}. #{query_msa.size} != #{ref_msa.size}")
				end
				#puts ">#{arr[2]}\n#{ref_msa}\n>#{arr[0]}\n#{query_msa}"
				puts "#{arr[2]}\t#{ref_msa}\t#{arr[0]}\t#{query_msa}\t#{cigar}\t#{arr[1]}"
				
				
			end
		end
	end

	def read_fasta(fasta : String)
		ref = {} of String => String
		id = ""
		seq = ""
		File.each_line(fasta) do |line|
			if line.match(/^>(\S+)/)
				unless id == ""
					ref[id] = seq.gsub(/\s/, "")
				end
				line =~ /^>(\S+)/
				id = $1
				seq = ""
				next
			else
				seq += line
			end
		end
		ref[id] = seq.gsub(/\s/, "")
		return ref
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
