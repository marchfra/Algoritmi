.PHONY: all
all:  ## Build all
	$(MAKE) -C ./ProjectileMotion


.PHONY: run
run:  ## Run all
	$(MAKE) run -C ./ProjectileMotion


.PHONY: plot
plot:  ## Create plot
	$(MAKE) plot -C ./ProjectileMotion


.PHONY: help
help:  ## Show this message
	@grep -E '^[a-zA-Z0-9_-]+:.*?## .*$$' $(MAKEFILE_LIST) \
	| sed -n 's/^\(.*\): \(.*\)## \(.*\)/\1|||\3/p' \
	| column -t  -s '|||'
