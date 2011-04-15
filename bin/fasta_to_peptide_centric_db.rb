#!/usr/bin/ruby

require 'rubygems'
require 'optparse'
require 'ms/ident/peptide/db'

Ms::Ident::Peptide::Db.cmdline(ARGV)
