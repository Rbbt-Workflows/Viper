require 'rbbt-util'
require 'rbbt/workflow'
require 'rbbt/sources/organism'

Misc.add_libdir if __FILE__ == $0

#require 'rbbt/sources/Viper'

module Viper
  extend Workflow

  input :modules, :text, "TF modules to use as regulon"
  task :regulon => :tsv do |tf_modules|
    tsv = TSV.setup({}, :key_field => "Regulator", :fields => ["Target", "Mode", "Likelihood"], :type => :double)

    require 'rbbt/sources/organism'
    organism = "Hsa/feb2014"

    found = Set.new
    tf_modules = Open.read(tf_modules) if Misc.is_filename?(tf_modules)
    tf_modules = tf_modules.to_s(false)  if TSV === tf_modules
    TSV.traverse StringIO.new(tf_modules), :type => :array, :into => tsv do |line|
      next if line =~ /^#/
      tf, desc, *targets_raw = line.split("\t")
      targets = []
      signs = []

      next if found.include? tf
      found << tf

      targets_raw.each do |target|
        if m = target.match(/(.*)\[(.*)\]/)
          #target = gene2ens[m[1]]
          #next if target.nil?
          target = m[1]
          targets << target
          signs << m[2]
        else
          #target = gene2ens[target]
          #next if target.nil?
          targets << target
          signs << 1
        end
      end
      
      signs = nil if signs.uniq == [1]
      [tf, [targets, signs, [] ]]
    end
    tsv
  end

  dep :regulon
  input :data, :tsv, "Expression data"
  input :min_genes, :integer, "Minumum number of genes per TF", 10
  task :viper => :tsv do |data, min_genes|
    regulon = step(:regulon).load

    data = TSV.open data, :unnamed => true unless TSV === data

    Open.write(file('regulon'), regulon.to_s)

    begin
      data_key = data.key_field
      organism = data.namespace || "Hsa/feb2014"
      reg_key, count = Organism.guess_id organism, regulon.column("Target").values.flatten
      if reg_key and data_key != reg_key
        log :translate, "Translating: #{[data_key, reg_key] * " => "}" unless data_key == reg_key
        data = data.change_key(reg_key, :identifiers => Organism.identifiers(organism)) 
      end
    rescue
      Log.warn "Could not normalize data identifiers: #{$!.message}"
    end

    data = data.to_list{|v| Misc.mean(v.compact.collect{|v| v.to_f})}  if data.type == :double

    require 'rbbt/util/R'
    script =<<-EOF
rbbt.require('viper')

regulon.tsv = rbbt.tsv('#{file('regulon')}')

regulon = list()
min.genes = #{min_genes}
for (tf in rownames(regulon.tsv)){

  info = regulon.tsv[tf,]

  targets.str = info$Target

  targets = strsplit(targets.str, '\\\\|') [[1]]

  if (length(targets) > min.genes){

    likelihood.str = info$Likelihood
    if (is.null(likelihood.str) || is.na(likelihood.str) || likelihood.str == ""){ likelihood = rep(1, length(targets))}
    else{ likelihood = as.numeric(strsplit(likelihood.str, '\\\\|')[[1]]) }

    tfmode.str = info$Mode
    if (is.null(tfmode.str) || is.na(tfmode.str) || tfmode.str == ""){ tfmode = rep(1, length(targets))}
    else{ tfmode = as.numeric(strsplit(tfmode.str, '\\\\|')[[1]]) }
    names(tfmode) = targets

    regulon[[tf]] = list(tfmode=tfmode, likelihood=likelihood)
  }
}

data = viper(data, regulon)
    EOF
    Open.write(file('script'), script)

    begin
      data.R script
    rescue
      raise RbbtException, "Viper failed: #{$!.message}"
    end
  end

  dep :viper
  task :profile => :tsv do
    step(:viper).load.transpose("Sample")
  end

  dep :regulon
  input :data, :tsv, "Expression data"
  input :main, :array, "Main samples"
  input :contrast, :array, "Contrast samples"
  input :organism, :string, "Organism code", Organism.default_code("Hsa")
  input :min_genes, :integer, "Minumum number of genes per TF", 10
  task :msviper => :tsv do |data,main,contrast,organism, min_genes|
    regulon = step(:regulon).load

    raise ParameterException, "No samples in the main category" if main.nil? or main.empty?
    raise ParameterException, "Only on sample in the main category" if main.length == 1

    all_genes = regulon.values.flatten.uniq

    Open.write(file('regulon'), regulon.to_s)

    organism = Organism.default_code("Hsa") if organism.nil?

    reg_key, count = Organism.guess_id organism, regulon.keys
    data = TSV.open(data, :cast => :to_f) unless TSV === data
    data_key = data.key_field
    data = data.to_list{|v| Misc.mean v}

    if (not reg_key.nil?) and data_key != reg_key
      data.identifiers = Organism.identifiers(organism)
      data = data.change_key(reg_key)
    end

    data = data.select(all_genes)

    if contrast.nil?
      contrast = data.fields - main
    end

    total = data.fields.length
    min = [contrast.length, main.length].min

    require 'distribution'

    permutations = Math.factorial(total) / (Math.factorial(min) * Math.factorial(total - min))

    raise ParameterException, "Not enough permutations to compute statistics" if permutations < 100
    permutations -= 1

    permutations = 1000 if permutations > 1000

    require 'rbbt/util/R'
    script =<<-EOF
rbbt.require('viper')

regulon.tsv = rbbt.tsv('#{file('regulon')}')

main.samples = #{R.ruby2R main}
contrast.samples = #{R.ruby2R contrast}

main = as.matrix(data[,main.samples])
contrast = as.matrix(data[,contrast.samples])

min.genes = #{min_genes}
regulon = list()
for (tf in rownames(regulon.tsv)){

  info = regulon.tsv[tf,]

  targets.str = info$Target

  targets = strsplit(targets.str, '\\\\|') [[1]]

  if (length(targets) >= min.genes){
    likelihood.str = info$Likelihood
    if (is.null(likelihood.str) || is.na(likelihood.str) || likelihood.str == ""){ likelihood = rep(1, length(targets))}
    else{ likelihood = as.numeric(strsplit(likelihood.str, '\\\\|')[[1]]) }

    tfmode.str = info$Mode
    if (is.null(tfmode.str) || is.na(tfmode.str) || tfmode.str == ""){ tfmode = rep(1, length(targets))}
    else{ tfmode = as.numeric(strsplit(tfmode.str, '\\\\|')[[1]]) }
    names(tfmode) = targets

    regulon[[tf]] = list(tfmode=tfmode, likelihood=likelihood)
  }

}

signature <- rowTtest(main, contrast)
nullModel <- ttestNull(main, contrast, per=#{permutations})

sig = signature$statistic

res = msviper(sig, regulon, nullModel, minsize = min.genes, verbose=FALSE)

data = data.frame(p.value=res$es$p.value, NES=res$es$nes)

data
    EOF
    Open.write(file('script'), script)

    begin
      data = data.R script

      data.key_field = "Associated Gene Name"
      data.namespace = organism

      data
    rescue
      raise RbbtException, "Viper failed: #{$!.message}"
    end
  end

end

#require 'Viper/tasks/basic.rb'

#require 'rbbt/knowledge_base/Viper'
#require 'rbbt/entity/Viper'

