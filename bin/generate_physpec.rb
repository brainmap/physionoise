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
require 'rubygems'
require 'logger'
require 'optparse'
require 'yaml'
require 'metamri'


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


def files(opts)
  defaults = { :in => '.', :matching => /.*/, :only_largest => nil }
  opts = defaults.merge(opts)
  
  files = Dir.entries(opts[:in]).select { |f| f =~ opts[:matching] }
  return files if opts[:only_largest].nil?
  files.sort_by { |f| File.size(File.join(opts[:in], f)) }.last(opts[:only_largest]).sort_by { |f| f }
end


# Assumptions: Respiratory files are large and in order.
def epis_and_associated_phys_files(raw_data_dir, phys_data_dir)
  epi_pattern = /fMRI/i
  
  phys_file_patterns = {
    :cardiac_signal      => /^ECGData_epiRT/,
    :cardiac_trigger     => /^TrigEcg_epiRT/,
    :respiration_signal  => /^RespData_epiRT/,
    :respiration_trigger => /^TrigResp_epiRT/
  }
  
  # scan raw data, pull out epi runs, sorted by timestamp
  visit = VisitRawDataDirectory.new(raw_data_dir)
  visit.scan
  visit.datasets.each { |ds| $LOG.debug "#{ds.directory} -- #{ds.series_description}" }
  epi_runs = visit.datasets.select { |ds| ds.series_description =~ epi_pattern }.sort_by { |ds| ds.timestamp }
  $LOG.debug epi_runs
  
  nruns = epi_runs.size
  
  # find nruns largest phys files in each phys file group
  physfile_groups = Hash.new
  phys_file_patterns.map do |type, pattern|
    found_files = files :in => phys_data_dir, :matching => pattern, :only_largest => nruns
    physfile_groups[type] = found_files 
  end
  physfile_groups.each { |type,files| $LOG.debug "#{type}: #{files}" }
  
  #
  yamlout = epi_runs.map { Hash.new }
  nruns.times do |i|
    yamlout[i][:series_description] = epi_runs[i].series_description
    yamlout[i][:run_directory] = epi_runs[i].directory
    yamlout[i][:phys_directory] = phys_data_dir
    yamlout[i][:bold_reps] = epi_runs[i].raw_image_files.first.bold_reps
    yamlout[i][:rep_time] = epi_runs[i].raw_image_files.first.rep_time / 1000.0
    physfile_groups.each_pair do |type,files|
      yamlout[i][type] = files[i]
    end
  end
  $LOG.debug "YAML output:\n#{YAML::dump(yamlout)}"
  
  return yamlout
end



if __FILE__ == $0
  opts = parse_options
  matches = epis_and_associated_phys_files(opts[:raw], opts[:phys])
  puts YAML::dump(matches) if opts[:preview]
  unless opts[:yaml].nil?
    File.open(opts[:yaml], 'w' ) do |out|
      YAML.dump(matches, out)
    end
  end
end


