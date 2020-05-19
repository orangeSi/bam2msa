def read_fasta(fasta : String)
		ref = {} of String => String
		id = ""
		seq = ""
		File.each_line(fasta) do |line|
			if line.match(/^>/)
				raise "error: not support line: #{line} in #{fasta}\n" if line.match(/^>\s/)	
				unless id == ""
					raise "error: id #{id} occur more than one times in #{fasta} \n" if ref.has_key?(id)
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
		raise "error: id #{id} occur more than one times in #{fasta} \n" if ref.has_key?(id)
		ref[id] = seq.gsub(/\s/, "")
		return ref
end

