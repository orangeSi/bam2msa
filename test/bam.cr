ifile = ARGV[0]
puts "input is #{ifile}"
Process.run("samtools view #{ifile}", shell: true) do |proc|
	while line = proc.output.gets
		puts line
	end
	puts "here you get end the file"
end
