def read_fasta(fasta : String, chrs : Array = [] of String)
  ref = {} of String => String
  id = ""
  seq = ""
  chrs_size = chrs.size
  skip_chr = false
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
      if chrs_size >= 1
        skip_chr = (chrs.includes?(id)) ? false : true
      end
      id = "" if skip_chr
    elsif skip_chr == false
      seq += line
    end
  end
  unless id == ""
    raise "error: id #{id} occur more than one times in #{fasta} \n" if ref.has_key?(id)
    ref[id] = seq.gsub(/\s/, "")
  end
  raise "sorry, occur a bug: not get all #{chrs} for #{fasta}. If you give a github issue, I will appreciate it~\n" if chrs_size >= 1 && ref.keys.size != chrs_size
  return ref
end
