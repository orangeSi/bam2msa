def read_line_from_path_or_io(file_path : String | IO::FileDescriptor) # path or STDIN
  if file_path.is_a?(String)
    File.each_line(file_path) do |line|
      yield line
    end
  else
    file_path.each_line do |line|
      yield line
    end
  end
end

def read_fasta(fasta : String | IO::FileDescriptor, chrs : Array = [] of String)
  # puts "chrs=#{chrs}"
  # puts "#start read fasta #{fasta}"
  ref = {} of String => String
  id = ""
  seq = Array(String).new(1000)
  chrs_size = chrs.size
  skip_chr = false
  got_chrs_size = 0

  read_line_from_path_or_io(fasta) do |line|
    if line.starts_with?(">")
      #raise "error: not support line: #{line} in #{fasta}\n" if line.match(/^>\s/)
      unless id == ""
        raise "error: id #{id} occur more than one times in #{fasta} \n" if ref.has_key?(id)
        ref[id] = seq.join("")
      end
      id = line[1..].split()[0]
      seq = Array(String).new(1000)
      # puts "reading #{id}"
      if chrs_size >= 1
        skip_chr = (chrs.includes?(id)) ? false : true
      end
      if skip_chr
        id = ""
        break if got_chrs_size == chrs_size
      else
        got_chrs_size += 1
      end
    elsif skip_chr == false
      seq << line
    end
  end
  unless id == ""
    raise "error: id #{id} occur more than one times in #{fasta} \n" if ref.has_key?(id)
    ref[id] = seq.join("")
  end
  raise "sorry, occur a bug: not get all #{chrs} for #{fasta}. If you give a github issue, I will appreciate it~\n" if chrs_size >= 1 && ref.keys.size != chrs_size
  # puts "#end read fasta #{fasta}"
  return ref
end

def test
  p! ARGV[0]
  read_fasta(ARGV[0]).each do |name, seq|
    puts "#{seq.size}\t#{name}"
  end
end

#test()
