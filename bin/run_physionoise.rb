#!/usr/bin/env ruby
#
# == Author
#   Kristopher J. Kosmatka, kk4@medicine.wisc.edu
#
# == Copyright
#   Copyright (c) 2010 WADRC Imaging Core.
#
require 'logger'
require 'optparse'
require 'yaml'


def run_physionoise_from_cli(runs, ARGV)
  run_physionoise_on(runs, ARGV)
end


if __FILE__ == $0
  spec = ARGV.shift
  runs = YAML.load_file(spec)
  run_physionoise_from_cli(runs, ARGV)
end
