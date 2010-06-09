require 'yaml'
require 'metamri'

class Physionoise
  
  attr_reader :raw_data_dir, :phys_data_dir, :spec
  attr_accessor :config
  
  def initialize(raw_data_dir, phys_data_dir, config = {})
    @raw_data_dir = raw_data_dir
    @phys_data_dir = phys_data_dir

    default_config = {}
    default_config[:phys_file_patterns] = {
      # :cardiac_signal      => /^ECGData_epiRT/,
      # :cardiac_trigger     => /^TrigEcg_epiRT/,
      # :respiration_signal  => /^RespData_epiRT/i,
      # :respiration_trigger => /^TrigResp_epiRT/i
      :cardiac_signal      => /^PPGData_epiRT/i,
      :cardiac_trigger     => /^PPGTrig_epiRT/i,
      :respiration_signal  => /^RespData_epiRT/i,
      :respiration_trigger => /^RespTrig_epiRT/i
    }
    default_config[:epi_pattern] = /fMRI/i
    
    @config = default_config.merge(config)
  end
  
  # Assumptions: Respiratory files are large and in order.
  def epis_and_associated_phys_files
    return @spec unless @spec.nil?
    
    epi_runs = select_epis
    @nruns = epi_runs.size

    physfile_groups = select_phys_files
    
    @spec = build_spec(epi_runs, physfile_groups)
    
    $LOG.debug "YAML output:\n#{YAML::dump(@spec)}"

    return @spec
  end
  
  def select_epis
    # scan raw data, pull out epi runs, sorted by timestamp
    visit = VisitRawDataDirectory.new(@raw_data_dir)
    visit.scan
    visit.datasets.each { |ds| $LOG.debug "#{ds.directory} -- #{ds.series_description}" }
    epi_runs = visit.datasets.select { |ds| ds.series_description =~ @config[:epi_pattern] }.sort_by { |ds| ds.timestamp }
    $LOG.debug epi_runs
    return epi_runs
  end
  
  def select_phys_files
    # find nruns largest phys files in each phys file group
    physfile_groups = Hash.new
    @config[:phys_file_patterns].map do |type, pattern|
      found_files = files :in => @phys_data_dir, :matching => pattern, :only_largest => @nruns
      physfile_groups[type] = found_files 
    end
    physfile_groups.each { |type,files| $LOG.debug "#{type}: #{files}" }
    return physfile_groups
  end
  
  def files(opts)
    defaults = { :in => '.', :matching => /.*/, :only_largest => nil }
    opts = defaults.merge(opts)

    files = Dir.entries(opts[:in]).select { |f| f =~ opts[:matching] }
    return files if opts[:only_largest].nil?
    files.sort_by { |f| File.size(File.join(opts[:in], f)) }.last(opts[:only_largest]).sort_by { |f| f }
  end
  
  def build_spec(epi_runs, physfile_groups)
    spec = epi_runs.map { Hash.new }
    @nruns.times do |i|
      spec[i][:series_description] = epi_runs[i].series_description
      spec[i][:run_directory] = epi_runs[i].directory
      spec[i][:phys_directory] = @phys_data_dir
      spec[i][:bold_reps] = epi_runs[i].raw_image_files.first.bold_reps
      spec[i][:rep_time] = epi_runs[i].raw_image_files.first.rep_time / 1000.0
      physfile_groups.each_pair do |type,files|
        spec[i][type] = files[i]
      end
    end
    return spec
  end
  
  def to_yaml(filename)
    File.open(filename, 'w' ) do |f| 
      YAML.dump(epis_and_associated_phys_files, f)
    end
  end
  
  def self.run_physionoise_on(runs, opts = Array.new)
    physionoise_cmd = File.join(File.dirname(__FILE__), '..', 'bin', 'physionoise.py')
    
    runs.each do |run|
      cardsig = File.join run[:phys_directory], run[:cardiac_signal]
      cardtrig = File.join run[:phys_directory], run[:cardiac_trigger]
      respsig = File.join run[:phys_directory], run[:respiration_signal]
      resptrig = File.join run[:phys_directory], run[:respiration_trigger]
      prefix = run[:series_description].gsub(/ /,'_')
      tr = run[:rep_time]
      num_tr = run[:bold_reps]

      cmdfmt = "python #{physionoise_cmd} -c %s -o %s -r %s -p %s --TR %s --numTR %s"
      cmd = cmdfmt % [cardsig, cardtrig, respsig, prefix, tr, num_tr]
      cmd = "#{cmd} #{opts.join(' ')}"
      puts cmd
      system(cmd)
    end
  end
  
end