# gw_base
An environment to write and test Dockerized scripts to process GW &amp; EM data

# Requirements:
Dave's homecooked, Docker images
Adding a .env file at the root to contain docker-compose settings:

VOL == local volume to which to mount this application

VOL_DB == local volume to which to mount the MySQL database

VOL_DB_CONFIG == local volume to mount the MySQL config (relative to ./db_configuration)

DB_PWD == root MySQL password
