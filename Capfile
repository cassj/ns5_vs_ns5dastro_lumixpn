###
# Do some basic differential expression analysis on the expression data.
#
# I AM HERE!

require 'catpaws'

#generic settings
set :aws_access_key,  ENV['AMAZON_ACCESS_KEY']
set :aws_secret_access_key , ENV['AMAZON_SECRET_ACCESS_KEY']
set :ec2_url, ENV['EC2_URL']
set :ssh_options, { :user => "ubuntu", :keys=>[ENV['EC2_KEYFILE']]}
set :key, ENV['EC2_KEY']
set :key_file, ENV['EC2_KEYFILE']
set :ami, 'ami-52794c26'  #EC2 eu-west-1 32bit Lucid
set :instance_type, 'm1.small'
set :s3cfg, ENV['S3CFG'] #location of ubuntu s3cfg file
set :working_dir, '/mnt/work'

set :nhosts, 1
set :group_name, 'ns5_vs_ns5dastro_lumixpn'

set :snap_id, `cat SNAPID`.chomp #ec2 eu-west-1 
set :vol_id, `cat VOLUMEID`.chomp #empty until you've created a new volume
set :ebs_size, 5  #give us enough space to generate some analysis data
set :ebs_zone, 'eu-west-1b'  #wherever the ami is
set :dev, '/dev/sdf'
set :mount_point, '/mnt/data'



#make a new EBS volume from this snap 
#cap EBS:create

#and mount your EBS
#cap EBS:attach
#cap EBS:format_xfs
#cap EBS:mount_xfs



set :xlsfile, 'raw_data_dot_CSV_ver3-same_as_finalreport.xls'


desc "Upload original xls file"
task :upload_data, :roles => group_name do
  upload(xlsfile, "#{mount_point}/#{xlsfile}")
end
before 'upload_data', 'EC2:start'


desc "create csv files from xls file"
task :make_csv, :roles => group_name do
  sudo "apt-get -y install libspreadsheet-parseexcel-perl"
  run "mkdir -p #{working_dir}/scripts"
  upload "scripts/convert_excel.pl", "#{working_dir}/scripts/convert_excel.pl"
  run "chmod +x #{working_dir}/scripts/convert_excel.pl"
  run "cd #{mount_point} && #{working_dir}/scripts/convert_excel.pl #{xlsfile}"
end
before 'make_csv', 'EC2:start'


#if you want to keep the results

#cap EBS:snapshot

#and then shut everything down:

# cap EBS:unmount
# cap EBS:detach
# cap EBS:delete - unless you're planning to use it again.
# cap EC2:stop




