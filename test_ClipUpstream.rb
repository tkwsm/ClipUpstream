#!/usr/bin/ruby

# require 'minitest/unit'
require 'minitest/autorun'
require './ClipUpstream.rb'
require '/home/wadaken/bin/GFF3Parser.rb'

class TC_ClipUpstream < MiniTest::Unit::TestCase

  class Test_ClipUpstream
    include ClipUpstream
  end

  def setup
    path2speAgff = "./spur_sample.gff"
    path2speAfas = "./spur_sample.fas"
    path2speBgff = "./hpul_sample.gff"
    path2speBfas = "./hpul_sample.fas"
    @speAgffh    = GFF3Parser.main_parser2( path2speAgff )
    @speAscagffh = GFF3Parser.create_scagffh( @speAgffh )
    @speAscah    = GFF3Parser.main_parser3( path2speAfas )
    @speBgffh    = GFF3Parser.main_parser2( path2speBgff )
    @speBscagffh = GFF3Parser.create_scagffh( @speBgffh )
    @speBscah    = GFF3Parser.main_parser3( path2speBfas )
    @speAscaid   = "Scaffold881"
    @speBscaid   = "scaffold1647"
    @cu = Test_ClipUpstream.new
    @strand          = "+"
    @scaffold_length = 100000
    @length_upstream = 10000
    @gstart          = 22000
    @gend            = 31000
    @region_start    = 12000
    @region_end      = 25000
    @max_gene_length = 3000
    @sample_genes    = [ "SPU_015325gn" ]
  end

  def test_cals_start
    tmpresult = @cu.calc_start( @strand, @scaffold_length, @length_upstream, 
                                @gstart, @gend, @max_gene_length )
    assert_equal( [ @region_start, @region_end ], tmpresult )
  end

  def test_get_genes_in_region
    tmpresult = @cu.get_genes_in_region( @speAscagffh, @speAgffh, @speAscaid, 
                                         @region_start, @region_end )
    assert_equal( @sample_genes, tmpresult )
  end

  def test_sampletest
    assert_equal( "good result", @cu.sampletest)
  end

end
