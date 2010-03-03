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


def run_physionoise_on(runs, opts)
  runs.each do |run|
    cardsig = File.join run[:phys_directory], run[:cardiac_signal]
    cardtrig = File.join run[:phys_directory], run[:cardiac_trigger]
    respsig = File.join run[:phys_directory], run[:respiration_signal]
    resptrig = File.join run[:phys_directory], run[:respiration_trigger]
    prefix = run[:series_description].gsub(/ /,'_')
    tr = run[:rep_time]
    num_tr = run[:bold_reps]
    
    cmdfmt = "python physionoise.py -c %s -o %s -r %s -p %s --TR %s --numTR %s"
    cmd = cmdfmt % [cardsig, cardtrig, respsig, prefix, tr, num_tr]
    cmd = "#{cmd} #{opts.join(' ')}"
    puts cmd
    system(cmd)
  end
end


if __FILE__ == $0
  spec = ARGV.shift
  runs = YAML.load_file(spec)
  run_physionoise_on(runs, ARGV)
end