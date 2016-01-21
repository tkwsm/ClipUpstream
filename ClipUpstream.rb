#!/usr/bin/ruby
#
 require '/home/wadaken/bin/GFF3Parser.rb'
# require './GFF3Parser.rb'

## listf format. you can use "-" as Gene name for an annonimous gene.
# 1  Wnt1   HpuGene1    SpuGene23
# 2  Wnt2   HpuGene21   SpuGene3999
# 3  -      HpuGene31   SpuGene12
# ....
#
if ARGV.size == 0
  print "USAGE: script.rb <listf> <speAgff> <speAfas> <speBgff> <speAfas>\n"
  exit
end

listf = open( ARGV.shift )
path2speAgff = ARGV.shift
path2speAfas = ARGV.shift
path2speBgff = ARGV.shift
path2speBfas = ARGV.shift
length_upst  = ( ARGV.shift ).to_i
max_gene_length  = ( ARGV.shift ).to_i
outdir       = ARGV.shift

list_h = {}
speAgid = ""
speBgid = ""
a = []
listf.each do |x|
  a = x.chomp.split("\s")
  listid  = a[0]
  gname   = a[1]
  speAgid = a[2]
  speAgid = speAgid.slice(/(\S+)-mRNA/, 1) if speAgid =~ /-mRNA/
  speBgid = a[3]
  speBgid = speBgid.slice(/(\S+)-mRNA/, 1) if speBgid =~ /-mRNA/
  list_h[ "#{listid}_#{gname}"] = [speAgid, speBgid]
end

############

def calc_start( strand, scaffold_length, length_upstream, gstart, gend, max_gene_length )
  region_start = 0
  region_end   = 0
  if    strand == "+"
    if gend - gstart <= max_gene_length
      region_end = gend
    else
      region_end = gstart + max_gene_length
    end
    if gstart - length_upstream > 0
      region_start = gstart - length_upstream
    else
      region_start = 1
    end
  elsif strand == "-"
    if gend - gstart < max_gene_length
      region_start = gstart
    else
      region_start = gend - max_gene_length
    end
    if scaffold_length > gend + length_upstream
      region_end = gend + length_upstream
    else
      region_end = scaffold_length
    end
  end
  return [ region_start, region_end ]
end

def get_genes_in_region( scagffh, gffh, scaid, region_start, region_end )
  genes  = []
  gstart = 0
  gend   = 0
  scagffh[scaid].each do |gid|
    gstart = gffh[ gid ]["gene"][3].to_i
    gend   = gffh[ gid ]["gene"][4].to_i
    if    gend < region_start
    elsif region_end < gstart
    else
      genes << gid
    end
  end
  return genes
end

def get_vista_annotation( general_order, gffh, genes, region_start, region_end )
  gorder = general_order
  annostr = ""
  annotation_lines_set = []
  annotation_lines = []
  tid = ""
  genes.each do |gid|
    str    = gffh[gid]["gene"][6]
    annostr = ">" if str == "+" and gorder == "+"
    annostr = "<" if str == "-" and gorder == "+"
    annostr = "<" if str == "+" and gorder == "-"
    annostr = ">" if str == "-" and gorder == "-"
    gffh[gid].each_key do |tid|
      annotation_lines = []
      next if tid == "gene"
# STDERR.puts "defline #{gffh[gid][tid]["mRNA"]}"
      if    gorder == "+"     
        tstart = gffh[gid][tid]["mRNA"][0][3].to_i - region_start
        tend   = gffh[gid][tid]["mRNA"][0][4].to_i - region_start
        tstart = 1 if tstart < 0
        annotation_lines << "#{annostr} #{tstart} #{tend} #{tid}\n"
        if gffh[gid][tid]["CDS"] == []
          gffh[gid][tid]["exon"].each do |exoninfo|
            ftype  = exoninfo[2]
            estart = exoninfo[3].to_i - region_start
            eend   = exoninfo[4].to_i - region_start
            next if estart < -1 or eend < -1
            annotation_lines << "#{estart} #{eend} #{ftype}\n"
          end     
        else
          gffh[gid][tid]["CDS"].each do |exoninfo|
            ftype  = "exon"
            estart = exoninfo[3].to_i - region_start
            eend   = exoninfo[4].to_i - region_start
            next if estart < -1 or eend < -1
            annotation_lines << "#{estart} #{eend} #{ftype}\n"
          end     
          gffh[gid][tid]["3-utr"].each do |exoninfo|
            ftype  = "utr"
            estart = exoninfo[3].to_i - region_start
            eend   = exoninfo[4].to_i - region_start
            next if estart < -1 or eend < -1
            annotation_lines << "#{estart} #{eend} #{ftype}\n"
          end     
          gffh[gid][tid]["5-utr"].each do |exoninfo|
            ftype  = "utr"
            estart = exoninfo[3].to_i - region_start
            eend   = exoninfo[4].to_i - region_start
            next if estart < -1 or eend < -1
            annotation_lines << "#{estart} #{eend} #{ftype}\n"
          end     
        end
      elsif gorder == "-"
        tstart = region_end - gffh[gid][tid]["mRNA"][0][4].to_i + 1
        tend   = region_end - gffh[gid][tid]["mRNA"][0][3].to_i + 1
        tstart = 1 if tstart < 0
        annotation_lines << "#{annostr} #{tstart} #{tend} #{tid}\n"
        if gffh[gid][tid]["CDS"] == []
          gffh[gid][tid]["exon"].each do |exoninfo|
            ftype  = exoninfo[2]
            estart = region_end - exoninfo[4].to_i + 1
            eend   = region_end - exoninfo[3].to_i + 1
            next if estart < -1 or eend < -1
#            annotation_lines << "#{estart} #{eend} #{ftype} #{exoninfo[3]} #{exoninfo[4]} #{region_start}\n"
            annotation_lines << "#{estart} #{eend} #{ftype}\n"
          end     
        else
          gffh[gid][tid]["CDS"].each do |exoninfo|
            ftype  = "exon"
            estart = region_end - exoninfo[4].to_i + 1
            eend   = region_end - exoninfo[3].to_i + 1
            next if estart < -1 or eend < -1
#            annotation_lines << "#{estart} #{eend} #{ftype} #{exoninfo[3]} #{exoninfo[4]} #{region_start}\n"
            annotation_lines << "#{estart} #{eend} #{ftype}\n"
          end     
          gffh[gid][tid]["3-utr"].each do |exoninfo|
            ftype  = "utr"
            estart = region_end - exoninfo[4].to_i + 1
            eend   = region_end - exoninfo[3].to_i + 1
            next if estart < -1 or eend < -1
#            annotation_lines << "#{estart} #{eend} #{ftype} #{exoninfo[3]} #{exoninfo[4]} #{region_start}\n"
            annotation_lines << "#{estart} #{eend} #{ftype}\n"
          end     
          gffh[gid][tid]["5-utr"].each do |exoninfo|
            ftype  = "utr"
            estart = region_end - exoninfo[4].to_i + 1
            eend   = region_end - exoninfo[3].to_i + 1
            next if estart < -1 or eend < -1
#            annotation_lines << "#{estart} #{eend} #{ftype} #{exoninfo[3]} #{exoninfo[4]} #{region_start}\n"
            annotation_lines << "#{estart} #{eend} #{ftype}\n"
          end     
        end
      end
      annotation_lines[1..-1] = annotation_lines[1..-1].sort{|x, y| x.split("\s")[0].to_i <=> y.split("\s")[0].to_i }
      annotation_lines << "\n"
      annotation_lines_set << annotation_lines
    end
  end
  annotation_lines_set.sort!{|x, y| x[0].split("\s")[1].to_i <=> y[0].split("\s")[1].to_i}
  return annotation_lines_set
end

############

speAgffh    = GFF3Parser.main_parser2( path2speAgff )
speAscagffh = GFF3Parser.create_scagffh( speAgffh )
speAscah    = GFF3Parser.main_parser3( path2speAfas )
speBgffh    = GFF3Parser.main_parser2( path2speBgff )
speBscagffh = GFF3Parser.create_scagffh( speBgffh )
speBscah    = GFF3Parser.main_parser3( path2speBfas )

list_h.each_key do |k|

  outprefix = k
  speAoutfa = File.new("#{outdir}/#{outprefix}.1.fas", "w+")
  speBoutfa = File.new("#{outdir}/#{outprefix}.2.fas", "w+")
  speAoutan = File.new("#{outdir}/#{outprefix}.1.ann", "w+")
  speBoutan = File.new("#{outdir}/#{outprefix}.2.ann", "w+")

  speAgid, speBgid = list_h[ k ]
STDERR.puts speAgid
STDERR.puts speAgffh[ speAgid ]
  speA_scID = speAgffh[ speAgid ]["gene"][0]
  speA_gstt = speAgffh[ speAgid ]["gene"][3].to_i
  speA_gend = speAgffh[ speAgid ]["gene"][4].to_i
  speA_gstr = speAgffh[ speAgid ]["gene"][6]

  speB_scID = speBgffh[ speBgid ]["gene"][0]
  speB_gstt = speBgffh[ speBgid ]["gene"][3].to_i
  speB_gend = speBgffh[ speBgid ]["gene"][4].to_i
  speB_gstr = speBgffh[ speBgid ]["gene"][6]

  speA_scln  = speAscah[ speA_scID ].length
  speB_scln  = speBscah[ speB_scID ].length
  speA_regst, speA_regen = calc_start(speA_gstr, speA_scln, length_upst, speA_gstt, speA_gend, max_gene_length )
  speB_regst, speB_regen = calc_start(speB_gstr, speB_scln, length_upst, speB_gstt, speB_gend, max_gene_length)

  speA_genes = get_genes_in_region( speAscagffh, speAgffh, speA_scID, speA_regst, speA_regen )
  speB_genes = get_genes_in_region( speBscagffh, speBgffh, speB_scID, speB_regst, speB_regen )

  if speA_gstr == "+"
    speAoutfa.print speAscah[ speA_scID ].subseq(speA_regst, speA_regen).to_fasta( "#{speA_scID}:plus:#{speA_regst}-#{speA_regen}", "60" )
  else
    speAoutfa.print speAscah[ speA_scID ].subseq(speA_regst, speA_regen).complement.to_fasta( "#{speA_scID}:minus:#{speA_regen}-#{speA_regst}", "60" )
  end

  if speB_gstr == "+"
    speBoutfa.print speBscah[ speB_scID ].subseq(speB_regst, speB_regen).to_fasta( "#{speB_scID}:plus:#{speB_regst}-#{speB_regen}", "60" )
  else
    speBoutfa.print speBscah[ speB_scID ].subseq(speB_regst, speB_regen).complement.to_fasta( "#{speB_scID}:minus:#{speB_regen}-#{speB_regst}", "60" )
  end

  speA_annotation_lines_set = get_vista_annotation( speA_gstr, speAgffh, speA_genes, speA_regst, speA_regen )
  speA_annotation_lines_set.each do |speA_annotation_lines| 
    speA_annotation_lines.each{ |line| speAoutan.print line }
  end
  speB_annotation_lines_set = get_vista_annotation( speB_gstr, speBgffh, speB_genes, speB_regst, speB_regen )
  speB_annotation_lines_set.each do |speB_annotation_lines| 
    speB_annotation_lines.each{ |line| speBoutan.print line }
  end

#  print "speA scaffoldID: #{speA_scID}\t"
#  print "speA gene start: #{speA_gstt}\t"
#  print "speA gene end  : #{speA_gend}\t"
#  print "speA gene strand: #{speA_gstr}\n\n"
#
#  print "speB scaffoldID: #{speB_scID}\t"
#  print "speB gene start: #{speB_gstt}\t"
#  print "speB gene end  : #{speB_gend}\t"
#  print "speB gene strand: #{speB_gstr}\n\n"
#
#  print "speA scaffold-length: #{ speA_scln }\n"
#  print "speA gene region: #{speA_regst} #{speA_regen}\n"
#  print "speB scaffold-length: #{ speB_scln }\n"
#  print "speB gene region: #{speB_regst} #{speB_regen}\n"

end

