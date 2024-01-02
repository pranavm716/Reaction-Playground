<template>
  <div :id="jsmeContainerName"></div>
  <input :id="smilesContainerName" hidden="hidden" v-model="smiles">
</template>

<script>
import { defineComponent } from 'vue'
export default defineComponent(
  {
    name: 'ChemDrawTool',
    props: ['identifier'],
    emits: ['smiles'],
    data () {
      return {
        smiles: ''
      }
    },
    computed: {
      jsmeContainerName () {
        return 'jsmeContainer-' + this.identifier
      },
      smilesContainerName () {
        return 'smilesContainer-' + this.identifier
      },
      storeSmilesFunctionName () {
        return 'storeSmiles_' + this.identifier
      }
    },
    mounted () {
      if (document.getElementById('my-jsme')) return // was already loaded
      const jsmeScript = document.createElement('script')
      jsmeScript.src = 'jsme/jsme.nocache.js'
      jsmeScript.type = 'text/javascript'
      jsmeScript.language = 'javascript'
      jsmeScript.id = 'my-jsme'
      document.head.appendChild(jsmeScript)

      const jsmeOnLoad = document.createElement('script')
      jsmeOnLoad.innerHTML = `
        function jsmeOnLoad() {
            jsmeApplet = new JSApplet.JSME("${this.jsmeContainerName}", "380px", "340px");
            jsmeApplet.setCallBack("AfterStructureModified", ${this.storeSmilesFunctionName});
        }
      `
      document.head.appendChild(jsmeOnLoad)

      const storeSmiles = document.createElement('script')
      storeSmiles.innerHTML = `
        function ${this.storeSmilesFunctionName}(event) {
          smiles = event.src.smiles();
          const smilesContainer = document.getElementById("${this.smilesContainerName}");
          smilesContainer.value = smiles;
          smilesContainer.dispatchEvent(new Event('input'));
        }
      `
      document.head.appendChild(storeSmiles)
    },
    watch: {
      smiles (value) {
        this.smiles = value
        this.$emit('smiles', this.smiles)
      }
    }
  }
)
</script>
