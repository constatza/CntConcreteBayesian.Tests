/*****************************************************************************
   Copyright 2018 The TensorFlow.NET Authors. All Rights Reserved.

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
******************************************************************************/

using Google.Protobuf;
using System;

namespace Tensorflow
{
    internal sealed class SessionOptions : IDisposable
    {
        public SafeSessionOptionsHandle Handle { get; }

        public SessionOptions(string target = "", ConfigProto config = null)
        {
            Handle = c_api.TF_NewSessionOptions();
            c_api.TF_SetTarget(Handle, target);
            if (config != null)
                SetConfig(config);
        }

        public void Dispose()
            => Handle.Dispose();

        private unsafe void SetConfig(ConfigProto config)
        {
            var bytes = config.ToByteArray();

            fixed (byte* proto2 = bytes)
            {
                using (var status = new Status())
                {
                    c_api.TF_SetConfig(Handle, (IntPtr)proto2, (ulong)bytes.Length, status.Handle);
                    status.Check(false);
                }
            }
        }
    }
}
