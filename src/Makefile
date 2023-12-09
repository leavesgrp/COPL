CC       := gcc
CFLAGS   := -O

BuildDir := ../build
Target   := copl

Source = $(wildcard *.c)
OBJS   = $(Source:%.c=$(BuildDir)/%.o)

all: $(BuildDir) | $(Target)

$(Target): $(OBJS)
	@$(CC) -o $(BuildDir)/$(Target) $(OBJS) -lm

$(BuildDir)/%.o: %.c
	@$(CC) $(CFLAGS) -c $< -o $@

clean:
	@rm -rf $(BuildDir)

$(BuildDir):
	@mkdir $@

