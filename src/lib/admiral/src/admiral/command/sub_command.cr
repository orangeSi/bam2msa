abstract class Admiral::Command
  private module Run
    def run
      raise Error.new "Invalid subcommand: #{@argv[0]?}"
    end
  end

  # Invokes a sub command by name, passing `self` as the parent.
  abstract def sub(command, *args, **params)

  private macro inherited
    macro finished
      \{% if !SubCommands::SPECS.empty? %}
        include Run
      \{% end %}
    end

    private struct SubCommands
      SPECS = {} of String => NamedTuple(
        type: Admiral::Command.class,
        description: Tuple(String, String?),
        short: String?
      )

      macro finished
        def locate
          \{% for name, spec in SPECS %}
            if \{{ name }} == @name || \{{ spec[:short] }} == @name
              return \{{ spec[:type].id }}
            end
          \{% end %}
        end
      end

      def self.locate(name)
        new(name.to_s).locate
      end

      def self.invoke(name, *args, **params)
        new(name).invoke(*args, **params)
      end

      def initialize(name : ::Admiral::StringValue)
        initialize name.value
      end

      def initialize(@name : String)
      end

      def invoke(*args, **params)
        if sub_command_class = locate
          sub_command_class.run(*args, **params, program_name: @name)
        else
          raise ::Admiral::Error.new("Invalid subcommand: #{@name}.")
        end
      end
    end

    def sub(command, *args, **params)
      SubCommands.invoke(command, *args, **params, parent: self)
    end
  end

  # Registers a subcommand.
  #
  # ```crystal
  # # hello.cr
  # class Hello < Admiral::Command
  #   class Planetary < Admiral::Command
  #     def run
  #       puts "Hello World"
  #     end
  #   end
  #
  #   class Municipality < Admiral::Command
  #     def run
  #       puts "Hello Denver"
  #     end
  #   end
  #
  #   register_subcommand planet : Planetary
  #   register_subcommand city : Municipality
  #
  #   def run
  #     puts help
  #   end
  # end
  #
  # HelloWorld.run
  # ```
  #
  # ```sh
  # $ crystal build ./hello.cr
  # $ ./hello planet
  # Hello World
  # $ ./hello city
  # Hello Denver
  # ```
  macro register_sub_command(command, type = nil, *, description = nil, short = nil)
    {%
      name = command.is_a?(TypeDeclaration) ? command.var.stringify : command.id.stringify
      type = command.is_a?(TypeDeclaration) ? command.type : type

      SubCommands::SPECS[name] = {
        type: type,
        description: {
          name + (short ? ", #{short.id}" : ""),
          description
        },
        short: short
      }
    %}
  end
end
