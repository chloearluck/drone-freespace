#!/usr/bin/env ruby

# usage: .vtk2mp3 robot.vtk obstacle.vtk
# creates ivmodels, robot file, and scene file 
# (you need to manually add to .conf file to set a start and end configuration)

$MPK_PATH = "/Users/chloearluck/git/mpk_v1.0"

rob_filename = ARGV[0]
obs_filename = ARGV[1]


def createIVmodel infilename
  infile = File.readlines(infilename)
  i = 4
  while (!infile[i].include? "POINTS")
    i = i+1
  end
  nverts = infile[i].split[1].to_i
  i = i+1

  verts = []
  while verts.length < nverts
    if infile[i].strip != ""
      verts.push(infile[i].strip)
    end
    i=i+1
  end

  while (!infile[i].include? "POLYGONS")
    i = i+1
  end
  nfaces = infile[i].split[1].to_i
  i = i+1

  faces = []
  while faces.length < nfaces
    if infile[i].strip != ""
      str = infile[i].strip
      str[0] = ''
      faces.push(str)
    end
    i=i+1
  end

  basename = File.basename(infilename, ".vtk")
  puts "Writing to file: "+$MPK_PATH+"/ivmodels/"+basename+".iv"
  ivoutfile = File.open($MPK_PATH+"/ivmodels/"+basename+".iv", "w+")
  ivoutfile.puts "#Inventor V2.0 ascii\n\n"
  ivoutfile.puts "DEF "+basename+" mpkObstacle {"
  ivoutfile.puts "  Coordinate3 {\n    point ["

  for i in 0..(nverts-1) do
    ivoutfile.print "      "+verts[i]
    if i < (nverts-1)
      ivoutfile.print ","
    end
    ivoutfile.puts ""
  end
  ivoutfile.puts "    ]\n  }\n\n"

  ivoutfile.puts "  DEF __triangulate__ IndexedFaceSet {\n      coordIndex ["
  for i in 0..(nfaces-1) do
    ivoutfile.print "        "
    arr = faces[i].split
    arr.each do |n|
      ivoutfile.print n+", "
    end
    if i < (nfaces-1)
      ivoutfile.puts "-1,"
    else
      ivoutfile.puts "-1"
    end
  end

  ivoutfile.puts "      ]\n    }\n}"

  if (verts.length != nverts || faces.length != nfaces)
    puts "ERROR: problem ready vtk file "+infile
  end
end

createIVmodel rob_filename
createIVmodel obs_filename

robtext = '#MPK v1.0 robot

param {
  3 { cyclic }
}

joint x_transl {
  Transl1 1 0 0  -10  10
  param 0
}
joint y_transl {
  parent x_transl
  Transl1 0 1 0  -10  10
  param 1
}
joint z_transl {
  parent y_transl
  Transl1 0 0 1  -10  10
  param 2
}
joint z_rot {
  parent z_transl
  Rot1 0 0 1  -3.142  3.142
  param 3 
  model0 "BASENAME.iv"
  tracePoint
}'

rob_basename = File.basename(rob_filename, ".vtk")
robtext.sub!('BASENAME', rob_basename)
puts "Writing to file: "+$MPK_PATH+"/robots/"+rob_basename+".rob"
roboutfile = File.open($MPK_PATH+"/robots/"+rob_basename+".rob", "w+")
roboutfile.puts robtext

### Create scene file
obs_basename = File.basename(obs_filename, ".vtk")
puts "Writing to file: "+$MPK_PATH+"/scenes/"+rob_basename+"_"+obs_basename+".iv"
sceneoutfile = File.open($MPK_PATH+"/scenes/"+rob_basename+"_"+obs_basename+".iv", "w+")

sceneoutfile.puts "#Inventor V2.0 ascii\n\n"
sceneoutfile.puts "DEF robot mpkRobot {"
sceneoutfile.puts "  fileName \"#{rob_basename}.rob\""
sceneoutfile.puts "}"

sceneoutfile.puts "Separator {"
sceneoutfile.puts "  Material { transparency 0.5 }"
sceneoutfile.puts "  BaseColor { rgb 1 0 0 }"
sceneoutfile.puts "  DEF #{obs_basename} mpkIncludeFile {"
sceneoutfile.puts "    name \"#{obs_basename}.iv\""
sceneoutfile.puts "  }"
sceneoutfile.puts "}"

