using System;
using System.Collections.Generic;

namespace Godot
{
    internal static class MarshalUtils
    {
        private static Dictionary<object, object> ArraysToDictionary(object[] keys, object[] values)
        {
            var ret = new Dictionary<object, object>();

            for (var i = 0; i < keys.Length; i++)
            {
                ret.Add(keys[i], values[i]);
            }

            return ret;
        }

        private static void DictionaryToArrays(Dictionary<object, object> from, out object[] keysTo, out object[] valuesTo)
        {
            var keys = from.Keys;
            keysTo = new object[keys.Count];
            keys.CopyTo(keysTo, 0);

            var values = from.Values;
            valuesTo = new object[values.Count];
            values.CopyTo(valuesTo, 0);
        }

        private static Type GetDictionaryType()
        {
            return typeof(Dictionary<object, object>);
        }
    }
}
