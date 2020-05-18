t = "6S16M2D138M3D136M5D230M6I79M9I"
#ddt =~ /\d+[MIDNSHP=X]/
puts split_cigar(t)

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
		else
			clen = "#{clen}#{e}"
		end
	end
	return cigars
end
