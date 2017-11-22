require 'open3'
require 'rake/clean'

PROBLEM_NAME = "Hokudai"
ROUND_ID = 16981
TESTER = "#{PROBLEM_NAME}.jar"
SEED = 3
TEST_TYPE = ["random", "perfect"][0]

CLEAN.include %w(data/* *.gcda *.gcov *.gcno *.png)

desc 'c++ file compile'
task :compile do
  sh("g++ -std=c++11 -D__NO_INLINE__ -W -Wall -Wno-sign-compare -O2 -o #{PROBLEM_NAME} #{PROBLEM_NAME}.cpp")
end

desc 'exec and view result'
task run: [:compile] do
  sh("java -jar ./#{TESTER} -vis -seed #{SEED} -exec './#{PROBLEM_NAME}'")
end

desc 'check single'
task one: [:compile] do
  if ENV["debug"]
    sh("time java -jar #{TESTER} -seed #{SEED} -debug -novis -exec './#{PROBLEM_NAME}'")
  else
    sh("time ./#{PROBLEM_NAME} < testcases/testcase_#{SEED}.in > result.txt")
    sh("./score_evaluator.out testcases/testcase_#{SEED}.in result.txt")
  end
end

desc 'check for windows'
task windows: [:compile] do
  sh("java -jar ./#{TESTER} -novis -seed #{SEED} -exec ./#{PROBLEM_NAME}.exe")
end

desc 'check out of memory'
task :debug do
  sh("g++ -std=c++11 -W -Wall -g -fsanitize=address -fno-omit-frame-pointer -Wno-sign-compare -O2 -o #{PROBLEM_NAME} #{PROBLEM_NAME}.cpp")
  sh("time java -jar #{TESTER} -seed #{SEED} -novis -exec './#{PROBLEM_NAME}'")
end

desc 'check how many called each function'
task :coverage do
  sh("g++ -W -Wall -Wno-sign-compare -o #{PROBLEM_NAME} --coverage #{PROBLEM_NAME}.cpp")
  sh("time java -jar #{TESTER} -seed #{SEED} -novis -exec './#{PROBLEM_NAME}'")
end

desc "example test"
task sample: [:compile] do
  run_test(1..10)
end

desc "production test"
task test: [:compile] do
  run_test(1..100)
end

desc "system test"
task final: [:compile] do
  run_test(1001..2000)
end

desc "system test production"
task production: [:compile] do
  run_test(1..20000)
end

desc "check select seed"
task seeds: [:compile] do
  run_test(File.readlines("seeds.txt").map(&:to_i))
end

def run_test(seeds)
  File.open("result.txt", "w") do |file|
    seeds.each do |seed|
      puts "seed = #{seed}"
      file.puts("----- !BEGIN! ------")
      file.puts("Seed = #{seed}")

      sh("./#{PROBLEM_NAME} < testcases/testcase_#{seed}.in > output")
      data = Open3.capture3("./score_evaluator.out testcases/testcase_#{seed}.in output")
      file.puts(data.select { |d| d.is_a?(String) }.flat_map { |d| d.split("\n") })
      file.puts("----- !END! ------")
    end
  end

  ruby "scripts/analyze.rb #{seeds.size}"
end

task default: :compile

