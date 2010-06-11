#!/usr/bin/env ruby
#
# == Synopsis
#
# == Examples
#
# == Usage
#
# == Options
#
# == Author
#
# == Copyright
#   Copyright (c) 2010 WADRC Imaging Core.
#
$LOAD_PATH.unshift(File.dirname(__FILE__) + '/../lib')

require 'rubygems'
require 'logger'
require 'optparse'
require 'yaml'
require 'metamri'
require 'physiospec'


$LOG = Logger.new($stdout)
$LOG.level = Logger::ERROR


def parse_options
  options = Hash.new
  parser = OptionParser.new do |opts|
    opts.banner = "Usage: physpec.rb --raw=RAW_DATA_DIR --phys=PHYS_DATA_DIR [options]"
    opts.on("-r", "--raw RAW_DATA_DIR", "Specify the raw data directory") do |v|
      abort "* Raw data directory not found\n#{parser.help}" unless File.directory?(v)
      options[:raw] = v
    end
    opts.on("-p", "--phys PHYS_DATA_DIR", "Specify the physiology data directory") do |v|
      abort "* Physiology data directory not found\n#{parser.help}" unless File.directory?(v)
      options[:phys] = v
    end
    opts.on("-o", "--output OUTPUT_FILE", "Save the spec to a file") do |v|
      abort "* Output file already exists\n#{parser.help}" if File.exist?(v)
      options[:yaml] = v
    end
    opts.on("-e", "--preview", "Preview the spec without saving") do |v|
      options[:preview] = v
    end
    opts.on("-v", "--[no-]verbose", "Be wordy") do |v|
      $LOG.level = Logger::INFO
    end
    opts.on("-d", "--[no-]debug", "Display Debug messages") do |v|
      $LOG.level = Logger::DEBUG
    end
  end
  parser.parse!
  if options[:raw].nil? or options[:phys].nil?
    raise OptionParser::MissingArgument, "\n#{parser.help}"
  end

  $LOG.debug "Options:"
  options.each_pair { |k,v| $LOG.debug("    %-20s %s" % [k,v]) }
  $LOG.debug "Arguments"
  ARGV.each { |arg| $LOG.debug("    %s" % arg) }
  
  return options
end


if __FILE__ == $0
  opts = parse_options
  physiospec = Physionoise.new(opts[:raw], opts[:phys])
  puts YAML::dump(physiospec.epis_and_associated_phys_files) if opts[:preview]
  unless opts[:yaml].nil?
    physiospec.to_yaml(opts[:yaml])
  end
end
