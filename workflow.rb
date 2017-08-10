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
    TSV.traverse StringIO.new(tf_modules), :type => :array, :into => tsv do |line|
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
  task :viper => :tsv do |data|
    regulon = step(:regulon).load

    data = TSV.open data unless TSV === data

    Open.write(file('regulon'), regulon.to_s)

    begin
      data_key = data.key_field
      organism = data.namespace || "Hsa/feb2014"
      reg_key, count = Organism.guess_id organism, regulon.column("Target").values.flatten
      log :translate, "Translating: #{[data_key, reg_key] * " => "}" unless data_key == reg_key
      data = data.change_key(reg_key, :identifiers => Organism.identifiers(organism)) if reg_key and data_key != reg_key
    rescue
      Log.warn "Could not normalize data identifiers: #{$!.message}"
    end

    data = data.to_list{|v| Misc.mean(v)}  if data.type == :double

    require 'rbbt/util/R'
    script =<<-EOF
library(viper)

regulon.tsv = rbbt.tsv('#{file('regulon')}')

regulon = list()
for (tf in rownames(regulon.tsv)){

  info = regulon.tsv[tf,]

  targets.str = info$Target

  targets = strsplit(targets.str, '\\\\|') [[1]]

  likelihood.str = info$Likelihood
  if (is.null(likelihood.str) || is.na(likelihood.str)){ likelihood = rep(1, length(targets))}
  else{ likelihood = as.numeric(strsplit(likelihood.str, '\\\\|')[[1]]) }

  tfmode.str = info$Mode
  if (is.null(tfmode.str) || is.na(tfmode.str)){ tfmode = rep(1, length(targets))}
  else{ tfmode = as.numeric(strsplit(tfmode.str, '\\\\|')[[1]]) }
  names(tfmode) = targets

  regulon[[tf]] = list(likelihood=likelihood, tfmode=tfmode)
}

data = viper(data, regulon)
    EOF
    Open.write(file('script'), script)

    data.R script
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
  task :msviper => :tsv do |data,main,contrast,organism|
    regulon = step(:regulon).load

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

    require 'rbbt/util/R'
    script =<<-EOF
library(viper)

regulon.tsv = rbbt.tsv('#{file('regulon')}')

main.samples = #{R.ruby2R main}
contrast.samples = #{R.ruby2R contrast}

main = as.matrix(data[,main.samples])
contrast = as.matrix(data[,contrast.samples])

min.genes = 25
regulon = list()
for (tf in rownames(regulon.tsv)){

  info = regulon.tsv[tf,]

  targets.str = info$Target

  targets = strsplit(targets.str, '\\\\|') [[1]]

  if (length(targets) > min.genes){
    likelihood.str = info$Likelihood
    if (is.null(likelihood.str) || is.na(likelihood.str)){ likelihood = rep(1, length(targets))}
    else{ likelihood = as.numeric(strsplit(likelihood.str, '\\\\|')[[1]]) }

    tfmode.str = info$Mode
    if (is.null(tfmode.str) || is.na(tfmode.str)){ tfmode = rep(1, length(targets))}
    else{ tfmode = as.numeric(strsplit(tfmode.str, '\\\\|')[[1]]) }
    names(tfmode) = targets

    regulon[[tf]] = list(tfmode=tfmode, likelihood=likelihood)
  }

}

signature <- rowTtest(main, contrast)
nullModel <- ttestNull(main, contrast, per=1000)

sig = signature$statistic

res = msviper(sig, regulon, nullModel, minsize = min.genes, verbose=FALSE)

data = data.frame(p.value=res$es$p.value, NES=res$es$nes)

data
    EOF
    Open.write(file('script'), script)

    data = data.R script

    data.key_field = "Associated Gene Name"
    data.namespace = organism

    data

  end

end

#require 'Viper/tasks/basic.rb'

#require 'rbbt/knowledge_base/Viper'
#require 'rbbt/entity/Viper'

