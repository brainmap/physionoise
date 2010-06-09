#!/usr/bin/env ruby
#
# == Author
#   Kristopher J. Kosmatka, kk4@medicine.wisc.edu
#
# == Copyright
#   Copyright (c) 2010 WADRC Imaging Core.
#
$LOAD_PATH.unshift(File.join(File.dirname(__FILE__), '..', 'lib'))
require 'logger'
require 'optparse'
require 'yaml'
require 'physionoise'


def run_physionoise_from_cli(runs, *args)
  Physionoise.run_physionoise_on(runs, args)
end


if __FILE__ == $0
  spec = ARGV.shift
  runs = YAML.load_file(spec)
  run_physionoise_from_cli(runs, ARGV)
end
