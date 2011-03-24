#!/usr/bin/ruby

require 'rubygems'
require 'optparse'
require 'ms/ident/peptidedb'

Ms::Ident::Peptidedb.cmdline(ARGV)
