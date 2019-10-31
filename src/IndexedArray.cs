using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Landis.Extension.GiplDamm
{
    public partial class GiplDamm
    {
        private class IndexedArray<T>
        {
            private readonly T[] _vals;

            public IndexedArray(int minIndex, int maxIndex)
            {
                MinIndex = minIndex;
                MaxIndex = maxIndex;
                _vals = new T[maxIndex - minIndex + 1];
            }

            public T this[int index] { get { return _vals[index - MinIndex]; } set { _vals[index - MinIndex] = value; } }

            public int MinIndex { get; }
            public int MaxIndex { get; }

            public int Count => MaxIndex - MinIndex + 1;
            public T Max => _vals.Max();
            public T Min => _vals.Min();

            public IndexedArray<T> Copy()
            {
                var g = new IndexedArray<T>(MinIndex, MaxIndex);
                for (var i = MinIndex; i <= MaxIndex; ++i)
                    g[i] = this[i];

                return g;
            }
        }
    }
}
