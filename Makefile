build:
	docker-compose build

run-volco:
	docker-compose run volco python volco.py --gcode=$(GCODE) --sim=$(SIM) --printer=$(PRINTER)

test:
	docker-compose run volco pytest


.PHONY: build run-volco test
