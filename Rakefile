
require 'rubygems'
require 'rake'

require 'jeweler'
Jeweler::Tasks.new do |gem|
  # gem is a Gem::Specification... see http://docs.rubygems.org/read/chapter/20 for more options
  gem.name = "ms-error_rate"
  gem.homepage = "http://github.com/jtprince/ms-error_rate"
  gem.license = "MIT"
  gem.summary = %Q{An mspire library for calculating error rates in MS/MS identifications (FDRs).}
  gem.description = %Q{aids for creating and calculating error rates using target-decoy searches and sample validation.}
  gem.email = "jtprince@gmail.com"
  gem.authors = ["John Prince"]
  # Include your dependencies below. Runtime dependencies are required when using your gem,
  # and development dependencies are only needed for development (ie running rake tasks, tests, etc)
  #  gem.add_runtime_dependency 'jabber4r', '> 0.1'
  #  gem.add_development_dependency 'rspec', '> 1.2.3'
  gem.rubyforge_project = 'mspire'
  gem.add_runtime_dependency("ms-core", ">= 0.0.2")
  gem.add_runtime_dependency("ms-ident", ">= 0.0.20")
  gem.add_development_dependency "spec-more", ">= 0"
  gem.add_development_dependency "jeweler", "~> 1.5.2"
  gem.add_development_dependency "rcov", ">= 0"
end
Jeweler::RubygemsDotOrgTasks.new

require 'rake/testtask'
Rake::TestTask.new(:spec) do |spec|
  spec.libs << 'lib' << 'spec'
  spec.pattern = 'spec/**/*_spec.rb'
  spec.verbose = true
end

#require 'rcov/rcovtask'
#Rcov::RcovTask.new do |spec|
#  spec.libs << 'spec'
#  spec.pattern = 'spec/**/*_spec.rb'
#  spec.verbose = true
#end

task :default => :spec

require 'rake/rdoctask'
Rake::RDocTask.new do |rdoc|
  version = File.exist?('VERSION') ? File.read('VERSION') : ""

  rdoc.rdoc_dir = 'rdoc'
  rdoc.title = "ms-error_rate #{version}"
  rdoc.rdoc_files.include('README*')
  rdoc.rdoc_files.include('lib/**/*.rb')
end
