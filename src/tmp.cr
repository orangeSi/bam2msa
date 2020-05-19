
t = RG.new(100, 300)
puts t.s
puts t.e
struct RG
		property s, e
		def initialize(@s : Int32, @e : Int32)
		end
end
