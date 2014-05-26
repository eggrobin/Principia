using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using principia.integrators;

namespace ConsoleApplication1 {
  class Program {
    static void Main(string[] args) {
      for (int i = 1; i <= 10; ++i) {
        Console.WriteLine(System.DateTime.Now);
        principia.integrators.SPRKTests.Run();
        Console.WriteLine(System.DateTime.Now);
        Console.ReadLine();
      }
    }
  }
}
