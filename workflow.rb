require 'rbbt-util'
require 'rbbt/workflow'

Misc.add_libdir if __FILE__ == $0

#require 'rbbt/sources/Viper'

Workflow.require_workflow "CSBC"
module Viper
  extend Workflow

  input :modules, :text, "TF modules to use as regulon"
  task :regulon => :tsv do |tf_modules|
    tsv = TSV.setup({}, :key_field => "Regulator", :fields => ["Target", "Mode", "Likelihood"], :type => :double)

    require 'rbbt/sources/organism'
    organism = "Hsa/feb2014"

    found = Set.new
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

    Open.write(file('regulon'), regulon.to_s)

    data_key = data.key_field
    reg_key, count = Organism.guess_id "Hsa/feb2014", regulon.keys
    data = data.change_key(reg_key).to_list{|v| Misc.mean(v)} if data_key != reg_key

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

str(regulon)

data = viper(data, regulon)
    EOF
    Open.write(file('script'), script)

    data.R script
  end

  dep :viper
  task :profile => :tsv do
    step(:viper).load.transpose("Sample")
  end

end

#require 'Viper/tasks/basic.rb'

#require 'rbbt/knowledge_base/Viper'
#require 'rbbt/entity/Viper'

