#!/usr/bin/env ruby
infilename = ARGV[0]
outfilename = ARGV[1]

infile = File.readlines(infilename)
outfile = File.open(outfilename, "w")

nverts = infile[0].split[0].to_i
nfaces = infile[0].split[1].to_i

outfile.puts("# vtk DataFile Version 3.0\n")
outfile.puts("vtk output\n")
outfile.puts("ASCII\n")
outfile.puts("DATASET POLYDATA\n")
outfile.puts("POINTS #{nverts} double\n")

i = 2
while (infile[i] != "\n")
  outfile.puts(infile[i])
  i=i+1
end

outfile.puts("\nPOLYGONS #{nfaces} #{nfaces*4}\n")

i=i+1
while (infile[i] != "\n")
  outfile.puts("#{infile[i].split.size} " + infile[i])
  i=i+1
end

outfile.close()
